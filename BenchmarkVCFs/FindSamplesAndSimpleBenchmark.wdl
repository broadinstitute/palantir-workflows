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
                reference_files=[paired_data.left.vcf],
                reference_indices=[paired_data.left.index],
                haplotype_map=haplotype_map,
                check_all_file_pairs=check_all_file_pairs,
                fail_on_mismatch=fail_on_mismatch,
                check_only_matching_sample_names=check_only_matching_sample_names,
                crosscheck_by=crosscheck_by,
                lod_threshold=lod_threshold
        }

        # If the two files were a fingerprint match, run Benchmark
        if (length(MatchFingerprints.all_matched_pairs) > 0) {
            call SimpleBenchmark.SimpleBenchmark as SimpleBenchmark {
                input:
                    base_vcf=paired_data.left.vcf,
                    base_vcf_index=paired_data.left.index,
                    base_output_sample_name=paired_data.left.output_name,
                    base_vcf_sample_name=paired_data.left.sample_name,
                    query_vcf=paired_data.right.vcf,
                    query_vcf_index=paired_data.right.index,
                    query_output_sample_name=paired_data.right.output_name,
                    query_vcf_sample_name=paired_data.right.sample_name,
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    evaluation_intervals=paired_data.left.eval_intervals,
                    score_field=score_field,
                    experiment=experiment,
                    extra_column_name=extra_column_name,
                    extra_column_value=extra_column_value,
                    create_igv_session=create_igv_sessions,
                    igv_session_name=igv_session_name,
                    preemptible=preemptible
            }
        }
    }

#    scatter(pair in MatchFingerprints.all_matched_pairs) {
#        Array[File] array_fingerprint_pairs = [pair.left, pair.right]
#    }
#
#    Array[File] flat_fingerprint_pairs = flatten(array_fingerprint_pairs)
#
#    scatter(paired_data in cross(base_vcf_data, query_vcf_data)) {
#        call MatchVcfData {
#            input:
#                fingerprint_matched_pairs=flat_fingerprint_pairs,
#                base_vcf_data=paired_data.left,
#                query_vcf_data=paired_data.right
#        }
#    }
#
#    scatter(pair_status in MatchVcfData.sample_pair_status) {
#        if (pair_status.right) {
#            Pair[BaseVcfData, VcfData] matched_vcf_data = pair_status.left
#        }
#    }
#
#    Array[Pair[BaseVcfData, VcfData]] all_matched_vcf_data = select_all(matched_vcf_data)
#
#    scatter(matched_samples in all_matched_vcf_data) {
#        call SimpleBenchmark.SimpleBenchmark as SimpleBenchmark {
#            input:
#                base_vcf=matched_samples.left.vcf,
#                base_vcf_index=matched_samples.left.index,
#                base_output_sample_name=matched_samples.left.output_name,
#                base_vcf_sample_name=matched_samples.left.sample_name,
#                query_vcf=matched_samples.right.vcf,
#                query_vcf_index=matched_samples.right.index,
#                query_output_sample_name=matched_samples.right.output_name,
#                query_vcf_sample_name=matched_samples.right.sample_name,
#                ref_fasta=ref_fasta,
#                ref_index=ref_index,
#                stratifier_intervals=stratifier_intervals,
#                stratifier_labels=stratifier_labels,
#                evaluation_intervals=matched_samples.left.eval_intervals,
#                score_field=score_field,
#                experiment=experiment,
#                extra_column_name=extra_column_name,
#                extra_column_value=extra_column_value,
#                create_igv_session=create_igv_sessions,
#                igv_session_name=igv_session_name,
#                preemptible=preemptible
#        }
#    }

    output {
        Array[File] fingerprint_files = flatten(MatchFingerprints.fingerprint_files)

        Array[File] benchmark_summaries = select_all(SimpleBenchmark.SimpleSummary)
        Array[File] indel_stats = select_all(SimpleBenchmark.IndelDistributionStats)
        Array[File] snp_stats = select_all(SimpleBenchmark.SNPSubstitutionStats)
        Array[File] roc_stats = select_all(SimpleBenchmark.ROCStats)

        Array[File] igv_sessions = select_all(SimpleBenchmark.igv_session)
    }
}

task MatchVcfData {
    input {
        Array[File] fingerprint_matched_pairs    # Assumes matched pairs adjacent in Array, as computed in main workflow block
        BaseVcfData base_vcf_data
        VcfData query_vcf_data
    }

    command <<<
        set -xueo pipefail

        python3 << CODE
        import filecmp
        import gzip

        matched_files = ["~{sep="\", \"" fingerprint_matched_pairs}"]

        def get_file_contents(path):
            if path.split('.')[-1] == '.gz':
                with gzip.open(path, 'r') as file:
                    return file.read()
            else:
                with open(path, 'r') as file:
                    return file.read()

        base_contents = get_file_contents("~{base_vcf_data.vcf}")
        query_contents = get_file_contents("~{query_vcf_data.vcf}")

        # Returns true iff base file matches i and query file matches j bytewise
        def compare(i, j):
            return base_contents == get_file_contents(i) and query_contents == get_file_contents(j)

        with open("results.txt", "w") as file:
            if any([compare(x[0], x[1]) for x in zip(matched_files[0::2], matched_files[1::2])]):
                file.write("true")
            else:
                file.write("false")

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + 200 + " HDD"
        cpu: 4
        memory: 8 + "GB"
    }

    output {
        Boolean is_match = read_boolean("results.txt")
        Pair[Pair[BaseVcfData, VcfData], Boolean] sample_pair_status = ((base_vcf_data, query_vcf_data), is_match)
    }
}