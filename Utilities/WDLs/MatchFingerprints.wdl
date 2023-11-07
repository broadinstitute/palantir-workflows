version 1.0

workflow MatchFingerprints {
    input {
        Array[File] input_files
        Array[File] reference_files

        File haplotype_map

        Boolean check_all_file_pairs = true
        Boolean fail_on_mismatch = false
        Boolean check_only_matching_sample_names = false

        String crosscheck_by = "FILE"    # Or: READGROUP, LIBRARY, SAMPLE
        Float lod_threshold = -5
    }

    if (check_all_file_pairs) {
        scatter (pair in cross(input_files, reference_files)) {
            call CheckFingerprints as CheckAllFingerprints {
                input:
                    input_file=pair.left,
                    reference_file=pair.right,
                    haplotype_map=haplotype_map,
                    fail_on_mismatch=fail_on_mismatch,
                    check_only_matching_sample_names=check_only_matching_sample_names,
                    crosscheck_by=crosscheck_by,
                    lod_threshold=lod_threshold
            }
        }
    }

    if (!check_all_file_pairs) {
        scatter (pair in zip(input_files, reference_files)) {
            call CheckFingerprints as CheckCorrespondingFingerprints {
                input:
                    input_file=pair.left,
                    reference_file=pair.right,
                    haplotype_map=haplotype_map,
                    fail_on_mismatch=fail_on_mismatch,
                    check_only_matching_sample_names=check_only_matching_sample_names,
                    crosscheck_by=crosscheck_by,
                    lod_threshold=lod_threshold
            }
        }
    }

    # Collect all the matched pairs detected with GATK
    scatter (matched_pair in select_first([CheckAllFingerprints.sample_pair, CheckCorrespondingFingerprints.sample_pair])) {
        if (matched_pair.right) {
            Pair[File, File] matched_pairs = matched_pair.left
        }
    }

    output {
        Array[File] fingerprint_files = select_first([CheckAllFingerprints.fingerprint_file, CheckCorrespondingFingerprints.fingerprint_file])
        Array[Pair[File, File]] matched_pairs = select_all(matched_pairs)
    }
}

task CheckFingerprints {
    input {
        File input_file
        File reference_file

        File haplotype_map
        Boolean fail_on_mismatch
        Boolean check_only_matching_sample_names

        String crosscheck_by
        Float lod_threshold

        String output_name = "output"
        String gatk_tag = "4.4.0.0"
    }

    String crosscheck_mode = if check_only_matching_sample_names then "CHECK_SAME_SAMPLE" else "CHECK_ALL_OTHERS"

    command <<<
        set -xueo pipefail

        # Allow "UNEXPECTED_MATCH" at this stage using exit code 0 arg; otherwise causes exit code 1
        gatk CrosscheckFingerprints \
            -I ~{input_file} \
            -SI ~{reference_file} \
            -H ~{haplotype_map} \
            -O "~{output_name}.txt" \
            -LOD ~{lod_threshold} \
            --EXIT_CODE_WHEN_MISMATCH 0 \
            --CROSSCHECK_MODE ~{crosscheck_mode} \
            --CROSSCHECK_BY ~{crosscheck_by}

        # Check if any of the comparisons performed are classified as a MATCH
        if [ $(grep -v '#' "~{output_name}.txt" | awk '{print $3}' | grep "_MATCH" | wc -l) -gt 0 ];
        then
            echo "true" > result.txt
        else
            if [ "~{fail_on_mismatch}" = true ];
            then
                exit 1
            else
                echo "false" > result.txt
            fi
        fi

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        preemptible: 2
        disks: "local-disk " + ceil(size(input_file, "GB") + size(reference_file, "GB") + 20) + " HDD"
        memory: 8 + " GB"
    }

    output {
        Boolean contains_match = read_boolean("result.txt")
        File fingerprint_file = "~{output_name}.txt"
        Pair[Pair[File, File], Boolean] sample_pair = ((input_file, reference_file), contains_match)
    }
}
