version 1.0

struct IndexedFile {
    File main_file
    File index_file
}

workflow MatchFingerprints {
    input {
        Array[File] input_files
        Array[File] input_indices
        Array[File] reference_files
        Array[File] reference_indices

        File haplotype_map

        Boolean check_all_file_pairs = true
        Boolean fail_on_mismatch = false
        Boolean check_only_matching_sample_names = false

        String crosscheck_by = "FILE"    # Or: READGROUP, LIBRARY, SAMPLE
        Float lod_threshold = -5
    }

    scatter (file_pair in zip(input_files, input_indices)) {
        IndexedFile indexed_input_files = {"main_file": file_pair.left, "index_file": file_pair.right}
    }

    scatter (file_pair in zip(reference_files, reference_indices)) {
        IndexedFile indexed_reference_files = {"main_file": file_pair.left, "index_file": file_pair.right}
    }

    if (check_all_file_pairs) {
        Array[Pair[IndexedFile, IndexedFile]] cross_pair_iterates = cross(indexed_input_files, indexed_reference_files)
    }

    if (!check_all_file_pairs) {
        Array[Pair[IndexedFile, IndexedFile]] zip_pair_iterates = zip(indexed_input_files, indexed_reference_files)
    }

    Array[Pair[IndexedFile, IndexedFile]] pair_iterates = select_first([cross_pair_iterates, zip_pair_iterates])

    scatter (pair in pair_iterates) {
        call CheckFingerprints {
            input:
                input_file=pair.left.main_file,
                input_index=pair.left.index_file,
                reference_file=pair.right.main_file,
                reference_index=pair.right.index_file,
                haplotype_map=haplotype_map,
                fail_on_mismatch=fail_on_mismatch,
                check_only_matching_sample_names=check_only_matching_sample_names,
                crosscheck_by=crosscheck_by,
                lod_threshold=lod_threshold
        }
    }

    # Collect all the matched pairs detected with GATK
    scatter (matched_pair in CheckFingerprints.sample_pair) {
        if (matched_pair.right) {
            Pair[File, File] matched_pairs = matched_pair.left
        }
    }

    output {
        Array[File] fingerprint_files = CheckFingerprints.fingerprint_file
        Array[Pair[File, File]] all_matched_pairs = select_all(matched_pairs)
    }
}

task CheckFingerprints {
    input {
        File input_file
        File input_index
        File reference_file
        File reference_index

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
    Int memory = 8

    command <<<
        set -xueo pipefail

        # Create input index maps to handle cases when index file is not adjacent to main files
        echo -e "~{input_file}\t~{input_index}" > input_index_map.tsv
        echo -e "~{reference_file}\t~{reference_index}" > reference_index_map.tsv

        # Allow "UNEXPECTED_MATCH" at this stage using exit code 0 arg; otherwise causes exit code 1
        gatk --java-options "-Xmx~{memory-1}g" CrosscheckFingerprints \
            -I ~{input_file} \
            --INPUT_INDEX_MAP input_index_map.tsv \
            -SI ~{reference_file} \
            --SECOND_INPUT_INDEX_MAP reference_index_map.tsv \
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
        memory: memory + " GB"
    }

    output {
        Boolean contains_match = read_boolean("result.txt")
        File fingerprint_file = "~{output_name}.txt"
        Pair[Pair[File, File], Boolean] sample_pair = ((input_file, reference_file), contains_match)
    }
}
