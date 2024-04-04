version 1.0

struct IndexedFile {
    File main_file
    File index_file
    File extra_file    # For custom use cases to ensure main files stay "connected" to this when subsetting to matches
}

# Workflow for fingerprinting input_files against set of second_input_files
# Outputs list of files with matching samples inside. See below for full description of output type
workflow MatchFingerprints {
    input {
        Array[File] input_files
        Array[File] input_indices
        Array[File]? input_extra_files
        Array[File] second_input_files
        Array[File] second_input_indices
        Array[File]? second_input_extra_files

        File haplotype_map

        Boolean check_all_file_pairs = true    # Check all possible combinations of input files to fingerprint against each other
        Boolean fail_on_mismatch = false    # If any file compared fails to find a fingerprint match, fail the workflow
        Boolean check_only_matching_sample_names = false    # Check fingerprints only against matching sample names

        String crosscheck_by = "FILE"    # Or: READGROUP, LIBRARY, SAMPLE
        Float lod_threshold = -5
    }

    # Default "extra_file" to another copy of index if nothing provided
    Array[File] final_input_extra_files = select_first([input_extra_files, input_indices])
    Array[File] final_second_input_extra_files = select_first([second_input_extra_files, second_input_indices])

    scatter (file_pair in zip(zip(input_files, input_indices), final_input_extra_files)) {
        IndexedFile indexed_input_files = {"main_file": file_pair.left.left, "index_file": file_pair.left.right, "extra_file": file_pair.right}
    }

    scatter (file_pair in zip(zip(second_input_files, second_input_indices), final_second_input_extra_files)) {
        IndexedFile indexed_second_input_files = {"main_file": file_pair.left.left, "index_file": file_pair.left.right, "extra_file": file_pair.right}
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
                indexed_input_file=pair.left,
                indexed_second_input_file=pair.right,
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
            Pair[Pair[IndexedFile, IndexedFile], Array[Array[String]]] matched_pairs_and_samples = matched_samples
        }
    }

    output {
        Array[File] fingerprint_files = CheckFingerprints.fingerprint_file

        # Final output: Array of Pairs ("matches") which contain
        # left: a Pair of indexed input files with matching sample in given second_input file
        # right: an Array of Array[String] (each should be length = 2) taking an input sample name to a matching second_input sample name
        Array[Pair[Pair[IndexedFile, IndexedFile], Array[Array[String]]]] all_matched_pairs_and_samples = select_all(matched_pairs_and_samples)
    }
}

task CheckFingerprints {
    input {
        IndexedFile indexed_input_file
        IndexedFile indexed_second_input_file

        File haplotype_map
        Boolean fail_on_mismatch
        Boolean check_only_matching_sample_names

        String crosscheck_by
        Float lod_threshold

        String output_name = "output"
        String gatk_tag = "4.5.0.0"
    }

    parameter_meta {
        indexed_input_file: {
            localization_optional: true
        }
        indexed_second_input_file: {
            localization_optional: true
        }
    }

    String crosscheck_mode = if check_only_matching_sample_names then "CHECK_SAME_SAMPLE" else "CHECK_ALL_OTHERS"
    Int memory = 8

    command <<<
        set -xueo pipefail

        ## If inputs are VCFs then localize them
        # This avoids a Picard bug and requester-pays which hold our truth data VCFs we use often with this workflow
        # See https://github.com/broadinstitute/picard/issues/1927 for details
        if [ $(basename -s .vcf.gz "~{indexed_input_file.main_file}") != $(basename "~{indexed_input_file.main_file}") ]; then
            # Case where input is VCF
            gsutil cp ~{indexed_input_file.main_file} first_input.vcf.gz
            gsutil cp ~{indexed_input_file.index_file} first_input.vcf.gz.tbi
            TOOL_INPUT="first_input.vcf.gz"
            TOOL_INPUT_INDEX="first_input.vcf.gz.tbi"
        else
            # Not VCF input
            TOOL_INPUT="~{indexed_input_file.main_file}"
            TOOL_INPUT_INDEX="~{indexed_input_file.index_file}"
        fi
        if [ $(basename -s .vcf.gz "~{indexed_second_input_file.main_file}") != $(basename "~{indexed_second_input_file.main_file}") ]; then
            gsutil cp ~{indexed_second_input_file.main_file} second_input.vcf.gz
            gsutil cp ~{indexed_second_input_file.index_file} second_input.vcf.gz.tbi
            TOOL_SECOND_INPUT="second_input.vcf.gz"
            TOOL_SECOND_INPUT_INDEX="second_input.vcf.gz.tbi"
        else
            TOOL_SECOND_INPUT="~{indexed_second_input_file.main_file}"
            TOOL_SECOND_INPUT_INDEX="~{indexed_second_input_file.index_file}"
        fi

        # Create input index maps to handle cases when index file is not adjacent to main files
        echo -e "${TOOL_INPUT}\t${TOOL_INPUT_INDEX}" > input_index_map.tsv
        echo -e "${TOOL_SECOND_INPUT}\t${TOOL_SECOND_INPUT_INDEX}" > second_input_index_map.tsv

        # Allow "UNEXPECTED_MATCH" at this stage using exit code 0 arg; otherwise causes exit code 1
        gatk --java-options "-Xmx~{memory-1}g" CrosscheckFingerprints \
            -I ${TOOL_INPUT} \
            --INPUT_INDEX_MAP input_index_map.tsv \
            -SI ${TOOL_SECOND_INPUT} \
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
        disks: "local-disk " + ceil(size(indexed_input_file.main_file, "GB") + size(indexed_second_input_file.main_file, "GB") + 50) + " HDD"
        memory: memory + " GB"
    }

    output {
        File fingerprint_file = "~{output_name}.txt"
        Array[Array[String]] matching_sample_pairs = read_tsv("matching_sample_pairs.tsv")   # Returns [[""]] if empty

        # Output tuple using Pairs since mixed types: pair of given files and list of matching samples in each
        Pair[Pair[IndexedFile, IndexedFile], Array[Array[String]]] sample_pairs = ((indexed_input_file, indexed_second_input_file), matching_sample_pairs)
    }
}
