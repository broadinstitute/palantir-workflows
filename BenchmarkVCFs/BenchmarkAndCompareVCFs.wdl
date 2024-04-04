version 1.0

import "BenchmarkVCFs.wdl" as BenchmarkVCFs
import "CompareBenchmarks.wdl" as CompareBenchmarks

workflow BenchmarkAndCompareVCFs {
    input {
        Array[String] sample_ids
        Array[String] configurations

        Boolean include_counts = false
        Boolean generate_gc_plots = false

        Array[String]? order_of_samples
        Array[String]? order_of_configurations
        Array[Int]? deltas

        Int? preemptible
        # Variables on a input sample level
        Array[File] base_vcfs
        Array[File] base_vcf_indices
        Array[String] base_output_sample_names

        Array[File] evaluation_intervals

        Array[File] query_vcfs
        Array[File] query_vcf_indices
        Array[String] query_output_sample_names

        # Reference information
        File ref_fasta
        File ref_index

        # Subsetting inputs using intervals
        Array[File] stratifier_intervals = []
        Array[String] stratifier_labels = []

        # Evaluation inputs
        String score_field = "GQ"

        # Columns to add to output files
        String experiment = ""
        Array[String] extra_column_names = []
        Array[String] extra_column_values = []
    }


    scatter(i in range(length(base_vcfs))) {
        call BenchmarkVCFs.Benchmark as Benchmark {
            input:
                base_vcf=base_vcfs[i],
                base_vcf_index=base_vcf_indices[i],
                base_output_sample_name=base_output_sample_names[i],
                query_vcf=query_vcfs[i],
                query_vcf_index=query_vcf_indices[i],
                query_output_sample_name=query_output_sample_names[i],
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                score_field=score_field,
                experiment=experiment,
                extra_column_names=extra_column_names,
                extra_column_values=extra_column_values

        }
    }

    call CompareBenchmarks.CompareBenchmarks as CompareBenchmarks {
        input:
            benchmark_summaries = Benchmark.summary,
            sample_ids = sample_ids,
            configurations = configurations,
            stratifiers = stratLabels,
            include_counts = include_counts,
            generate_gc_plots = generate_gc_plots,
            order_of_samples = order_of_samples,
            order_of_configurations = order_of_configurations,
            deltas = deltas,
            preemptible = preemptible
    }

    scatter(i in range(length(sample_ids))) {
        call RenameSummary {
            input:
                input_summary = Benchmark.summary[i],
                suffix = sample_ids[i] + "_" + configurations[i],
                preemptible = preemptible
        }
    }

    output {
        File comparison_csv = CompareBenchmarks.comparison_csv
        File raw_data = CompareBenchmarks.raw_data
        Array[File]? gc_plots = CompareBenchmarks.gc_plots
        Array[File] summaries = RenameSummary.summary
    }
}

task RenameSummary {
    input {
        File input_summary
        String suffix
        Int? preemptible
    }

    command {
        mv ~{input_summary} ./summary_~{suffix}.csv
    }

    runtime {
        docker: "ubuntu:20.04"
        preemptible: select_first([preemptible, 0])
        memory: "2 GB"
        disks: "local-disk 20 HDD"
    }

    output{
        File summary = "summary_~{suffix}.csv"
    }
}