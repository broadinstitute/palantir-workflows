version 1.0

import "../BenchmarkVCFs/SimpleBenchmark.wdl" as BenchmarkVCFs

struct VcfFile {
    String sample_id
    File file
    File index
    Int num    # placeholder for position in list for cross product
}

struct TruthVcf {
    String sample_id
    File file
    File index
    File confidence_intervals
}

workflow NewFunctionalEquivalence {
    input {
        Array[String] sample_id
        Array[String] dataset

        Array[String] confidence_intervals

        Array[String] tool1_vcf
        Array[String] tool1_vcf_index
        Array[String] tool2_vcf
        Array[String] tool2_vcf_index

        Array[String] truth_vcf
        Array[String] truth_vcf_index

        Array[File]? stratifier_intervals
        Array[String]? stratifier_labels

        File ref_fasta
        File ref_index
        File hap_map

        String tool1_label
        String tool2_label
        String? additional_label

        Boolean signed_difference = false

        String gatkTag = "4.0.11.0"
        Boolean passingOnly = true
        Boolean requireMatchingGenotypes = true
        String vcfScoreField = "QUAL"
        Int? threadsVcfEval = 2
        Int? preemptible = 3
    }

    scatter (i in range(length(tool1_vcf))) {
        VcfFile tool1_inputs = {"sample_id": sample_id[i], "file": tool1_vcf[i], "index": tool1_vcf_index[i], "num": i}
    }

    scatter (i in range(length(tool2_vcf))) {
        VcfFile tool2_inputs = {"sample_id": sample_id[i], "file": tool2_vcf[i], "index": tool2_vcf_index[i], "num": i}
    }

    scatter (i in range(length(truth_vcf))) {
        TruthVcf truth_inputs = {"sample_id": sample_id[i], "file": truth_vcf[i], "index": truth_vcf_index[i], "confidence_intervals": confidence_intervals[i]}
    }

    # Evaluate against the truth files
    scatter (paired_vcfs in zip(tool1_vcf, truth_vcf)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalVsTruthTool1 {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=paired_vcfs.right.confidence_intervals,
                experiment="EvalVsTruthTool1"
        }
    }

    scatter (paired_vcfs in zip(tool2_vcf, truth_vcf)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalVsTruthTool1 {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=paired_vcfs.right.confidence_intervals,
                experiment="EvalVsTruthTool2"
        }
    }

    # Evaluate across the two tools
    scatter (paired_vcfs in zip(tool1_vcf, tool2_vcf)) {
        call BenchmarkVCFs.SimpleBenchmark as EvalInterTool {
            input:
                base_vcf=paired_vcfs.right.file,
                base_vcf_index=paired_vcfs.right.index,
                base_output_sample_name=paired_vcfs.right.sample_id,
                base_vcf_sample_name=paired_vcfs.right.sample_id,
                query_vcf=paired_vcfs.left.file,
                query_vcf_index=paired_vcfs.left.index,
                query_output_sample_name=paired_vcfs.left.sample_id,
                query_vcf_sample_name=paired_vcfs.left.sample_id,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                stratifier_intervals=stratifier_intervals,
                stratifier_labels=stratifier_labels,
                evaluation_intervals=paired_vcfs.right.confidence_intervals,
                experiment="EvalInterTool"
        }
    }

    # Evaluate within the same tool all possible pairs for both tools
    scatter (index in cross(range(length(tool1_vcf)), range(length(tool1_vcf)))) {
        if (index.left < index.right) {    # Only check when first has index less than second in cross product so no repeats
            call BenchmarkVCFs.SimpleBenchmark as EvalIntraTool1 {
                input:
                    base_vcf=tool1_vcf[index.left].file,
                    base_vcf_index=tool1_vcf[index.left].index,
                    base_output_sample_name=tool1_vcf[index.left].sample_id,
                    base_vcf_sample_name=tool1_vcf[index.left].sample_id,
                    query_vcf=tool1_vcf[index.right].file,
                    query_vcf_index=tool1_vcf[index.right].index,
                    query_output_sample_name=tool1_vcf[index.right].sample_id,
                    query_vcf_sample_name=tool1_vcf[index.right].sample_id,
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    experiment="EvalIntraTool1"
            }
        }
    }

    scatter (index in cross(range(length(tool2_vcf)), range(length(tool2_vcf)))) {
        if (index.left < index.right) {    # Only check when first has index less than second in cross product so no repeats
            call BenchmarkVCFs.SimpleBenchmark as EvalIntraTool1 {
                input:
                    base_vcf=tool2_vcf[index.left].file,
                    base_vcf_index=tool2_vcf[index.left].index,
                    base_output_sample_name=tool2_vcf[index.left].sample_id,
                    base_vcf_sample_name=tool2_vcf[index.left].sample_id,
                    query_vcf=tool2_vcf[index.right].file,
                    query_vcf_index=tool2_vcf[index.right].index,
                    query_output_sample_name=tool2_vcf[index.right].sample_id,
                    query_vcf_sample_name=tool2_vcf[index.right].sample_id,
                    ref_fasta=ref_fasta,
                    ref_index=ref_index,
                    stratifier_intervals=stratifier_intervals,
                    stratifier_labels=stratifier_labels,
                    experiment="EvalIntraTool2"
            }
        }
    }



    output {
#        Array[File] fe_plots = FEEvaluation.fe_plots
#        Array[File] f1_plots = F1Evaluation.f1_plots
#        Array[File] roc_plots = roc_plots_3
#        File merged_fe_plots = MergeFE.plots
#        File merged_f1_plots = MergeF1.plots
#        File merged_roc_plots = MergeROC.plots
#        File fe_summary = FEEvaluation.fe_summary
#        File f1_summary = F1Evaluation.f1_summary
#        Int fe_status = fe_status_combined
#        File html_report = CreateHTMLReport.report
    }
}