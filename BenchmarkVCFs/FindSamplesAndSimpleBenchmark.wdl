version 1.0

import "SimpleBenchmark.wdl" as SimpleBenchmark
import "../Utilities/WDLs/MatchFingerprints.wdl" as MatchFingerprints

struct VcfData {
    File vcf
    File index
    String output_name
    String sample_name
}

struct BaseVcfData {
    VcfData data
    File eval_intervals
}

workflow FindSamplesAndSimpleBenchmark {
    input {
        Array[File] base_vcfs
        Array[File] base_vcf_indices
        Array[String] base_sample_output_names
        Array[String] base_vcf_sample_names

        Array[File] query_vcfs
        Array[File] query_vcf_indices
        Array[String] query_sample_output_names
        Array[String] query_vcf_sample_names

        # Reference information
        File ref_fasta
        File ref_index

        # Subsetting inputs using intervals
        Array[File] stratifier_intervals = []
        Array[String] stratifier_labels = []

        # Evaluation inputs
        Array[File] evaluation_intervals
        String score_field = "GQ"

        # Columns to add to output files
        String? experiment
        String? extra_column_name
        String? extra_column_value

        # Fingerprint arguments
        File haplotype_map
        Boolean check_all_file_pairs = true
        Boolean fail_on_mismatch = false
        Boolean check_only_matching_sample_names = false
        String crosscheck_by = "FILE"    # Or: READGROUP, LIBRARY, SAMPLE
        Float lod_threshold = -5

        Boolean create_igv_sessions = false
        # Array[File]? optional_igv_bams    # Add this later
        String igv_session_name = "igv_session"

        # Toggle for more tries on preemptible machines for potentially cheaper runs at the risk of longer runtime
        Int preemptible = 3
    }

    call MatchFingerprints {
        input:
            input_files=query_vcfs,
            input_indices=query_vcf_indices,
            reference_files=base_vcfs,
            reference_indices=base_vcf_indices,
            haplotype_map=haplotype_map,
            check_all_file_pairs=check_all_file_pairs,
            fail_on_mismatch=fail_on_mismatch,
            check_only_matching_sample_names=check_only_matching_sample_names,
            crosscheck_by=crosscheck_by,
            lod_threshold=lod_threshold
    }

    scatter(base_data in zip(zip(base_vcfs, base_vcf_indices), zip(base_sample_output_names, base_vcf_sample_names))) {
        VcfData base_vcf_data = {"vcf": base_data.left.left, "index": base_data.left.right, "output_name": base_data.right.left, "sample_name": base_data.right.right}
    }

    scatter(eval_base_data in zip(base_vcf_data, evaluation_intervals)) {
        BaseVcfData eval_base_data = {"data": eval_base_data.left, "eval_intervals": eval_base_data.right}
    }

    scatter(query_data in zip(zip(query_vcfs, query_vcf_indices), zip(query_sample_output_names, query_vcf_sample_names))) {
        VcfData query_vcf_data = {"vcf": query_data.left.left, "index": query_data.left.right, "output_name": query_data.right.left, "sample_name": query_data.right.right}
    }

    scatter(paired_data in cross(eval_base_data, query_vcf_data)) {
        call MatchVcfData {
            input:
                fingerprint_matched_pairs=MatchFingerprints.all_matched_pairs,
                base_vcf_data=paired_data.left,
                query_vcf_data=paired_data.right
        }
    }

    scatter(pair_status in MatchVcfData.sample_pair_status) {
        if (pair_status.right) {
            Pair[BaseVcfData, VcfData] matched_vcf_data = pair_status.left
        }
    }

    Array[Pair[VcfData, VcfData]] all_matched_vcf_data = select_all(matched_vcf_data)

    scatter(matched_samples in all_matched_vcf_data) {
        call SimpleBenchmark {
            input:
                base_vcf=matched_samples.left.data.vcf,
                base_vcf_index=matched_samples.left.data.index,
                base_output_sample_name=matched_samples.left.data.output_name,
                base_vcf_sample_name=matched_samples.left.data.sample_name,
                query_vcf=matched_samples.right.vcf,
                query_vcf_index=matched_samples.right.index,
                query_output_sample_name=matched_samples.right.output_name,
                query_vcf_sample_name=matched_samples.right.sample_name,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=matched_samples.left.eval_intervals,
                score_field=score_field,
                experiment=experiment,
                extra_column_name=extra_column_name,
                extra_column_value=extra_column_value,
                create_igv_session=create_igv_sessions,
                igv_session_name=igv_session_name,
                preemptible=preemptible
        }
    }

    output {
        Array[File] fingerprint_files = MatchFingerprints.fingerprint_files

        Array[File] benchmark_summaries = SimpleBenchmark.SimpleSummary
        Array[File] indel_stats = SimpleBenchmark.IndelDistributionStats
        Array[File] snp_stats = SimpleBenchmark.SNPSubstitutionStats
        Array[File] roc_stats = SimpleBenchmark.ROCStats

        Array[File] igv_sessions = SimpleBenchmark.igv_session
    }
}

task MatchVcfData {
    input {
        Array[Pair[File, File]] fingerprint_matched_pairs
        BaseVcfData base_vcf_data
        VcfData query_vcf_data
    }

    command <<<
        set -xueo pipefail

        python3 << CODE

        matched_files = ["~{sep="\", \"" fingerprint_matched_pairs}"]
        print("Here is matched files")
        print(matched_files)

        with open("results.txt", "w") as file:
            search_string = "(~{base_vcf_data.vcf}, ~{query_vcf_data.vcf})"
            print(f"Search string is: {search_string}")
            if search_string in matched_files:
                file.write("true")
            else:
                file.write("false")

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + 200 + " HDD"
        cpu: 2
        memory: 4 + "GB"
    }

    output {
        Boolean is_match = read_boolean("results.txt")
        Pair[Pair[BaseVcfData, VcfData], Boolean] sample_pair_status = ((base_vcf_data, query_vcf_data), is_match)
    }
}