version 1.0

import "SimpleBenchmark.wdl" as SimpleBenchmark
import "../Utilities/WDLs/MatchFingerprints.wdl" as MatchFingerprints
import "../Utilities/WDLs/CombineTables.wdl" as CombineTables

struct VcfData {
    File vcf
    File index
    String output_name
    String sample_name
}

struct BaseVcfData {
    File vcf
    File index
    String output_name
    String sample_name
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

    scatter(base_data in zip(zip(zip(base_vcfs, base_vcf_indices), zip(base_sample_output_names, base_vcf_sample_names)), evaluation_intervals)) {
        BaseVcfData base_vcf_data = {"vcf": base_data.left.left.left, "index": base_data.left.left.right, "output_name": base_data.left.right.left, "sample_name": base_data.left.right.right, "eval_intervals": base_data.right}
    }

    scatter(query_data in zip(zip(query_vcfs, query_vcf_indices), zip(query_sample_output_names, query_vcf_sample_names))) {
        VcfData query_vcf_data = {"vcf": query_data.left.left, "index": query_data.left.right, "output_name": query_data.right.left, "sample_name": query_data.right.right}
    }

    scatter(paired_data in cross(base_vcf_data, query_vcf_data)) {
        call MatchFingerprints.MatchFingerprints as MatchFingerprints {
            input:
                input_files=[paired_data.right.vcf],
                input_indices=[paired_data.right.index],
                second_input_files=[paired_data.left.vcf],
                second_input_indices=[paired_data.left.index],
                second_input_extra_files=[paired_data.left.eval_intervals],
                haplotype_map=haplotype_map,
                check_all_file_pairs=check_all_file_pairs,
                fail_on_mismatch=fail_on_mismatch,
                check_only_matching_sample_names=check_only_matching_sample_names,
                crosscheck_by=crosscheck_by,
                lod_threshold=lod_threshold
        }
    }

    # If the two files had a fingerprint match, run Benchmark over the matching sample name pairs
    scatter(matched_files in flatten(MatchFingerprints.all_matched_pairs_and_samples)) {
        scatter(matched_sample_pair in matched_files.right) {
            call SimpleBenchmark.SimpleBenchmark as SimpleBenchmark {
                input:
                    base_vcf=matched_files.left.right.main_file,
                    base_vcf_index=matched_files.left.right.index_file,
                    base_output_sample_name=matched_sample_pair[1],
                    base_vcf_sample_name=matched_sample_pair[1],
                    query_vcf=matched_files.left.left.main_file,
                    query_vcf_index=matched_files.left.left.index_file,
                    query_output_sample_name=matched_sample_pair[0],
                    query_vcf_sample_name=matched_sample_pair[0],
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    evaluation_intervals=matched_files.left.right.extra_file,
                    score_field=score_field,
                    experiment=experiment,
                    extra_column_name=extra_column_name,
                    extra_column_value=extra_column_value,
                    check_fingerprint=false,
                    create_igv_session=create_igv_sessions,
                    igv_session_name=igv_session_name,
                    preemptible=preemptible
            }
        }
    }

    call CombineTables.CombineTables as CombineBenchmarkSummaries {
        input:
            tables=select_all(flatten(SimpleBenchmark.SimpleSummary))
    }

    call CombineTables.CombineTables as CombineSnpStats {
        input:
            tables=select_all(flatten(SimpleBenchmark.SNPSubstitutionStats))
    }

    call CombineTables.CombineTables as CombineIndelStats {
        input:
            tables=select_all(flatten(SimpleBenchmark.IndelDistributionStats))
    }

    call CombineTables.CombineTables as CombineRocStats {
        input:
            tables=select_all(flatten(SimpleBenchmark.ROCStats))
    }

    output {
        Array[File] fingerprint_files = flatten(MatchFingerprints.fingerprint_files)

        File benchmark_summaries = CombineBenchmarkSummaries.combined_table
        File snp_stats = CombineSnpStats.combined_table
        File indel_stats = CombineIndelStats.combined_table
        File roc_stats = CombineRocStats.combined_table

        Array[File] igv_sessions = select_all(flatten(SimpleBenchmark.igv_session))
    }
}
