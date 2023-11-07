version 1.0

workflow MatchFingerprints {
    input {
        Array[File] input_files
        Array[File]? input_indices
        Array[File] reference_files
        Array[File]? reference_indices

        File haplotype_map

        Boolean check_all_file_pairs = true
        Boolean fail_on_mismatch = false
        Boolean check_only_matching_sample_names = false

        String crosscheck_by = "FILE"    # Or: READGROUP, LIBRARY, SAMPLE
        Float lod_threshold = -5
    }

    if (check_all_file_pairs) {
        if (defined(input_indices) && (defined(reference_indices))) {
            scatter (pair in cross(zip(input_files, select_first([input_indices])), zip(reference_files, select_first([reference_indices])))) {
                call CheckFingerprints as CheckAllIndexFingerprints {
                    input:
                        input_file=pair.left.left,
                        input_index=pair.left.right,
                        reference_file=pair.right.left,
                        reference_index=pair.right.right,
                        haplotype_map=haplotype_map,
                        fail_on_mismatch=fail_on_mismatch,
                        check_only_matching_sample_names=check_only_matching_sample_names,
                        crosscheck_by=crosscheck_by,
                        lod_threshold=lod_threshold
                }
            }
        }

        if (!defined(input_indices) || !defined(reference_indices)) {
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
    }

    if (!check_all_file_pairs) {
        if (defined(input_indices) && (defined(reference_indices))) {
            scatter (pair in cross(zip(input_files, select_first([input_indices])), zip(reference_files, select_first([reference_indices])))) {
                call CheckFingerprints as CheckCorrespondingIndexFingerprints {
                    input:
                        input_file=pair.left.left,
                        input_index=pair.left.right,
                        reference_file=pair.right.left,
                        reference_index=pair.right.right,
                        haplotype_map=haplotype_map,
                        fail_on_mismatch=fail_on_mismatch,
                        check_only_matching_sample_names=check_only_matching_sample_names,
                        crosscheck_by=crosscheck_by,
                        lod_threshold=lod_threshold
                }
            }
        }

        if (!defined(input_indices) || !defined(reference_indices)) {
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
    }

    # Collect all the matched pairs detected with GATK
    scatter (matched_pair in select_first([CheckAllIndexFingerprints.sample_pair, CheckAllFingerprints.sample_pair, CheckCorrespondingIndexFingerprints.sample_pair, CheckCorrespondingFingerprints.sample_pair])) {
        if (matched_pair.right) {
            Pair[File, File] matched_pairs = matched_pair.left
        }
    }

    output {
        Array[File] fingerprint_files = select_first([CheckAllIndexFingerprints.fingerprint_file, CheckAllFingerprints.fingerprint_file, CheckCorrespondingIndexFingerprints.fingerprint_file, CheckCorrespondingFingerprints.fingerprint_file])
        Array[Pair[File, File]] all_matched_pairs = select_all(matched_pairs)
    }
}

task CheckFingerprints {
    input {
        File input_file
        File? input_index
        File reference_file
        File? reference_index

        File haplotype_map
        Boolean fail_on_mismatch
        Boolean check_only_matching_sample_names

        String crosscheck_by
        Float lod_threshold

        String output_name = "output"
        String gatk_tag = "4.4.0.0"
    }

    parameter_meta {
        input_file: {
            localization_optional: true
        }
        input_index: {
            localization_optional: true
        }
        reference_file: {
            localization_optional: true
        }
        reference_index: {
            localization_optional: true
        }
    }

    String crosscheck_mode = if check_only_matching_sample_names then "CHECK_SAME_SAMPLE" else "CHECK_ALL_OTHERS"

    command <<<
        set -xueo pipefail

        # Create input index maps to handle cases when index file is not adjacent to main files
        echo -e "~{input_file}\t~{input_index}" > input_index_map.tsv
        echo -e "~{reference_file}\t~{reference_index}" > reference_index_map.tsv

        # Allow "UNEXPECTED_MATCH" at this stage using exit code 0 arg; otherwise causes exit code 1
        gatk CrosscheckFingerprints \
            -I ~{input_file} \
            ~{true="--INPUT_INDEX_MAP input_index_map.tsv" false="" defined(input_index)} \
            -SI ~{reference_file} \
            ~{true="--SECOND_INPUT_INDEX_MAP reference_index_map.tsv" false="" defined(reference_index)} \
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
