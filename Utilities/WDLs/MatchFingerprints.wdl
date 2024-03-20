version 1.0

struct IndexedFile {
    File main_file
    File index_file
}

workflow MatchFingerprints {
    input {
        Array[File] input_files
        Array[File] input_indices
        Array[File] second_input_files
        Array[File] second_input_indices

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

    scatter (file_pair in zip(second_input_files, second_input_indices)) {
        IndexedFile indexed_second_input_files = {"main_file": file_pair.left, "index_file": file_pair.right}
    }

    if (check_all_file_pairs) {
        Array[Pair[IndexedFile, IndexedFile]] cross_pair_iterates = cross(indexed_input_files, indexed_second_input_files)
    }

    if (!check_all_file_pairs) {
        Array[Pair[IndexedFile, IndexedFile]] zip_pair_iterates = zip(indexed_input_files, indexed_second_input_files)
    }

    Array[Pair[IndexedFile, IndexedFile]] pair_iterates = select_first([cross_pair_iterates, zip_pair_iterates])

    scatter (pair in pair_iterates) {
        call CheckFingerprints {
            input:
                input_file=pair.left.main_file,
                input_index=pair.left.index_file,
                second_input_file=pair.right.main_file,
                second_input_index=pair.right.index_file,
                haplotype_map=haplotype_map,
                fail_on_mismatch=fail_on_mismatch,
                check_only_matching_sample_names=check_only_matching_sample_names,
                crosscheck_by=crosscheck_by,
                lod_threshold=lod_threshold
        }
    }

    # Collect all the matched pairs detected with GATK
    scatter (matched_samples in CheckFingerprints.sample_pairs) {
        if (length(select_first(matched_samples.right)) > 1) {
            Pair[Pair[File, File], Array[Array[String]]] matched_pairs_and_samples = matched_samples
        }
    }

    output {
        Array[File] fingerprint_files = CheckFingerprints.fingerprint_file

        # Final output: Array of Pairs ("matches") which contain
        # left: a Pair of input file with matching sample in given second_input file
        # right: an Array of Array[String] (each should be length = 2) taking an input sample name to a matching second_input sample name
        Array[Pair[Pair[File, File], Array[Array[String]]]] all_matched_pairs_and_samples = select_all(matched_pairs_and_samples)
    }
}

task CheckFingerprints {
    input {
        File input_file
        File input_index
        File second_input_file
        File second_input_index

        File haplotype_map
        Boolean fail_on_mismatch
        Boolean check_only_matching_sample_names

        String crosscheck_by
        Float lod_threshold

        String output_name = "output"
        String gatk_tag = "4.5.0.0"
    }

    parameter_meta {
        input_file: {
            localization_optional: true
        }
        input_index: {
            localization_optional: true
        }
        second_input_file: {
            localization_optional: true
        }
        second_input_index: {
            localization_optional: true
        }
    }

    String crosscheck_mode = if check_only_matching_sample_names then "CHECK_SAME_SAMPLE" else "CHECK_ALL_OTHERS"
    Int memory = 8

    command <<<
        set -xueo pipefail

        # Create input index maps to handle cases when index file is not adjacent to main files
        echo -e "~{input_file}\t~{input_index}" > input_index_map.tsv
        echo -e "~{second_input_file}\t~{second_input_index}" > second_input_index_map.tsv

        # Allow "UNEXPECTED_MATCH" at this stage using exit code 0 arg; otherwise causes exit code 1
        gatk --java-options "-Xmx~{memory-1}g" CrosscheckFingerprints \
            -I ~{input_file} \
            --INPUT_INDEX_MAP input_index_map.tsv \
            -SI ~{second_input_file} \
            --SECOND_INPUT_INDEX_MAP second_input_index_map.tsv \
            -H ~{haplotype_map} \
            -O "~{output_name}.txt" \
            -LOD ~{lod_threshold} \
            --EXIT_CODE_WHEN_MISMATCH 0 \
            --CROSSCHECK_MODE ~{crosscheck_mode} \
            --CROSSCHECK_BY ~{crosscheck_by}

        # Analyze table outputs and find matching sample pairs
        python << CODE
        import pandas as pd

        fp_df = pd.read_csv("~{output_name}.txt", sep="\t", comment="#")
        matching_rows = fp_df["RESULT"].apply(lambda x: "_MATCH" in x)
        fp_df[matching_rows][["LEFT_SAMPLE", "RIGHT_SAMPLE"]].to_csv("matching_sample_pairs.tsv", sep="\t", index=False, header=False)

        CODE

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        preemptible: 2
        disks: "local-disk " + ceil(size(input_file, "GB") + size(second_input_file, "GB") + 20) + " HDD"
        memory: memory + " GB"
    }

    output {
        File fingerprint_file = "~{output_name}.txt"
        Array[Array[String]] matching_sample_pairs = read_tsv("matching_sample_pairs.tsv")   # Returns [[""]] if empty

        # Output tuple using Pairs since mixed types: pair of given files and list of matching samples in each
        Pair[Pair[File, File], Array[Array[String]]] sample_pairs = ((input_file, second_input_file), matching_sample_pairs)
    }
}
