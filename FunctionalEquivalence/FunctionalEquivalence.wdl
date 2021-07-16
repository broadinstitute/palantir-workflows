version 1.0

import "../BenchmarkVCFs/BenchmarkVCFs.wdl" as BenchmarkVCFs
import "subworkflows/FEEvaluation.wdl" as FEEvaluation
import "subworkflows/F1Evaluation.wdl" as F1Evaluation
import "subworkflows/PlotROC.wdl" as PlotROC

workflow FunctionalEquivalence {
    input{
        Array[String] sample_id
        Array[String] dataset
        Array[String] replicate_no

        Array[String] confidence_intervals

        Array[String] tool1_vcf
        Array[String] tool1_vcf_index
        Array[String] tool2_vcf
        Array[String] tool2_vcf_index

        Array[String] truth_vcf
        Array[String] truth_vcf_index

        Array[File]? stratIntervals
        Array[String]? stratLabels

        File reference
        String referenceVersion = "hg38"
        File refIndex
        File refDict
        File hapMap

        String tool1_label
        String tool2_label
        String? additional_label

        Boolean signed_difference = false

        String gatkTag="4.0.11.0"
        Boolean passingOnly=true
        Boolean requireMatchingGenotypes=true
        String vcfScoreField = "QUAL"
        Int? threadsVcfEval=2
        Int? preemptible = 3
    }
        
    scatter (i in range(length(sample_id))) {
        if (truth_vcf[i] != "null") {
            call BenchmarkVCFs.Benchmark as EvalVsTruthTool1 {
                input:
                    confidenceInterval = confidence_intervals[i],
                    truthLabel = "truth",
                    evalLabel = dataset[i] + "." + replicate_no[i] + ".tool1",
                    truthVcf = truth_vcf[i],
                    truthVcfIndex = truth_vcf_index[i],
                    evalVcf = tool1_vcf[i],
                    evalVcfIndex = tool1_vcf_index[i],
                    stratIntervals = stratIntervals,
                    stratLabels = stratLabels,
                    reference = reference,
                    refIndex = refIndex,
                    refDict = refDict,
                    referenceVersion = referenceVersion,
                    hapMap = hapMap,
                    gatkTag = gatkTag,
                    passingOnly = passingOnly,
                    requireMatchingGenotypes = requireMatchingGenotypes,
                    vcfScoreField = vcfScoreField,
                    threadsVcfEval = threadsVcfEval,
                    preemptible = preemptible,
                    doIndelLengthStratification = false,
                    enableRefOverlap = true
            }

            call BenchmarkVCFs.Benchmark as EvalVsTruthTool2 {
                input:
                    confidenceInterval = confidence_intervals[i],
                    truthLabel = "truth",
                    evalLabel = dataset[i] + "." + replicate_no[i] + ".tool2",
                    truthVcf = truth_vcf[i],
                    truthVcfIndex = truth_vcf_index[i],
                    evalVcf = tool2_vcf[i],
                    evalVcfIndex = tool2_vcf_index[i],
                    stratIntervals = stratIntervals,
                    stratLabels = stratLabels,
                    reference = reference,
                    refIndex = refIndex,
                    refDict = refDict,
                    referenceVersion = referenceVersion,
                    hapMap = hapMap,
                    gatkTag = gatkTag,
                    passingOnly = passingOnly,
                    requireMatchingGenotypes = requireMatchingGenotypes,
                    vcfScoreField = vcfScoreField,
                    threadsVcfEval = threadsVcfEval,
                    preemptible = preemptible,
                    doIndelLengthStratification = false,
                    enableRefOverlap = true
            }

            Array[File] roc_tables_1 = flatten([
                select_all(EvalVsTruthTool1.snpRocs),
                select_all(EvalVsTruthTool1.nonSnpRocs),
                select_all(EvalVsTruthTool2.snpRocs),
                select_all(EvalVsTruthTool2.nonSnpRocs)
            ])
        
            call PlotROC.PlotROC as ROCPlot {
                input:
                    dataset_name = dataset[i],
                    files = roc_tables_1,
                    preemptible = preemptible
            }
            Array[File] roc_plots_1 = [ROCPlot.snp_plot, ROCPlot.snp_plot_zoomed, ROCPlot.snp_plot_zoomed_LCR, ROCPlot.indel_plot, ROCPlot.indel_plot_zoomed, ROCPlot.indel_plot_zoomed_LCR]
        } # if truth_vcf[i] != "null"

        call BenchmarkVCFs.Benchmark as VcfEval_Inter {
            input:
                evalLabel = "tool1" + dataset[i] + replicate_no[i],
                truthLabel = "tool2" + dataset[i] + replicate_no[i],
                evalVcf = tool1_vcf[i],
                evalVcfIndex = tool1_vcf_index[i],
                truthVcf = tool2_vcf[i],
                truthVcfIndex = tool2_vcf_index[i],
                confidenceInterval = confidence_intervals[i],
                reference = reference,
                refIndex = refIndex,
                refDict = refDict,
                hapMap = hapMap,
                stratIntervals = stratIntervals,
                stratLabels = stratLabels,
                referenceVersion = referenceVersion,
                threadsVcfEval = threadsVcfEval,
                gatkTag = gatkTag,
                requireMatchingGenotypes = requireMatchingGenotypes,
                passingOnly = passingOnly,
                vcfScoreField = vcfScoreField,
                preemptible = preemptible,
                doIndelLengthStratification = false,
                enableRefOverlap = true
        }

        call RenameSummary as RenameSummaryInter {
            input:
                input_summary = select_first([VcfEval_Inter.summary]),
                name = "inter" + dataset[i] + replicate_no[i],
                preemptible = preemptible
        }

        # Ideally we would want index j in [i+1;N] to avoid duplicates, but WDL does not support range to
        # start at values != 0. Therefore, start at zero, count until N - i - 1 and then always add i+1.
        #     j  in [i+1;N]     |  define j' = j-i-1
        # <=> j' in [0;N-i-1]
        scatter (j_minus_i_minus_1 in range(length(sample_id) - i - 1)) {
            Int j = j_minus_i_minus_1 + i + 1

            if (dataset[i] == dataset[j]) {
                call BenchmarkVCFs.Benchmark as VcfEval_IntraTool1 {
                    input:
                        evalLabel = "tool1" + dataset[i] + replicate_no[i],
                        truthLabel = "tool1" + dataset[j] + replicate_no[j],
                        evalVcf = tool1_vcf[i],
                        evalVcfIndex = tool1_vcf_index[i],
                        truthVcf = tool1_vcf[j],
                        truthVcfIndex = tool1_vcf_index[j],
                        confidenceInterval = confidence_intervals[i],
                        reference = reference,
                        refIndex = refIndex,
                        refDict = refDict,
                        hapMap = hapMap,
                        stratIntervals = stratIntervals,
                        stratLabels = stratLabels,
                        referenceVersion = referenceVersion,
                        threadsVcfEval = threadsVcfEval,
                        gatkTag = gatkTag,
                        requireMatchingGenotypes = requireMatchingGenotypes,
                        passingOnly = passingOnly,
                        vcfScoreField = vcfScoreField,
                        preemptible = preemptible,
                        doIndelLengthStratification = false,
                        enableRefOverlap = true
                }

                call RenameSummary as RenameSummaryIntraTool1 {
                    input:
                        input_summary = select_first([VcfEval_IntraTool1.summary]),
                        name = "tool1" + dataset[i] + replicate_no[i] + replicate_no[j],
                        preemptible = preemptible
                }

                call BenchmarkVCFs.Benchmark as VcfEval_IntraTool2 {
                    input:
                        evalLabel = "tool2" + dataset[i] + replicate_no[i],
                        truthLabel = "tool2" + dataset[j] + replicate_no[j],
                        evalVcf = tool2_vcf[i],
                        evalVcfIndex = tool2_vcf_index[i],
                        truthVcf = tool2_vcf[j],
                        truthVcfIndex = tool2_vcf_index[j],
                        confidenceInterval = confidence_intervals[i],
                        reference = reference,
                        refIndex = refIndex,
                        refDict = refDict,
                        hapMap = hapMap,
                        stratIntervals = stratIntervals,
                        stratLabels = stratLabels,
                        referenceVersion = referenceVersion,
                        threadsVcfEval = threadsVcfEval,
                        gatkTag = gatkTag,
                        requireMatchingGenotypes = requireMatchingGenotypes,
                        passingOnly = passingOnly,
                        vcfScoreField = vcfScoreField,
                        preemptible = preemptible,
                        doIndelLengthStratification = false,
                        enableRefOverlap = true
                }

                call RenameSummary as RenameSummaryIntraTool2 {
                    input:
                        input_summary = select_first([VcfEval_IntraTool2.summary]),
                        name = "tool2" + dataset[i] + replicate_no[i] + replicate_no[j],
                        preemptible = preemptible
                }
            } # if dataset[i] == dataset[j]

            # This is the horrible part of combining all inputs into one array.

            # Combine output from IntraTool1 (File?) and IntraTool2 (File?) to an Array[File] with 2 entries.
            Array[File] summaries1 = select_all([RenameSummaryIntraTool1.summary, RenameSummaryIntraTool2.summary])
        } # scatter j

        # summaries1 is now an Array[Array[File]] where the inner arrays consist of each either none or two Files
        # and the outer array originates from the scatter.
        Array[File] summaries2 = flatten(summaries1)

        # We also need to do this with the ROC tables. These are only generated if a truth_vcf[i] is present,
        # therefore roc_tables_1 is an Array[File]?. Either select that or return an empty array.
        Array[File] roc_tables_2 = select_first([roc_tables_1, []])

        # Very similar with the ROC plots. These are only generated if a truth_vcf[i] is present, therefore
        # roc_tables_1 is an Array[File]?. Either select that or return an empty array.
        Array[File] roc_plots_2 = select_first([roc_plots_1, []])
    } # scatter i

    # Analogous to assignment summaries2 = flatten(...). Here we combine all intra comparisons for samples[0...N]
    Array[File] summaries3 = flatten(summaries2)

    Array[File] inter_summaries = RenameSummaryInter.summary

    # At this point we want to add the intra summaries (summaries3) to the inter summaries (inter_summaries),
    # both of which are an Array[File]. Therefore, we define an array consisting of these two arrays and flatten it, thereby appending
    # one array to the other. This way, we have all files in one single array.
    Array[File] summaries4 = flatten([summaries3, inter_summaries])

    # We also need to flatten the roc_tables_2 array which has been created over the scatter i.
    Array[File] roc_tables_3 = flatten(roc_tables_2)

    call FEEvaluation.FEEvaluation {
        input:
            tool1_label = tool1_label,
            tool2_label = tool2_label,
            additional_label = additional_label,
            summaries = summaries4,
            preemptible = preemptible
    }

    call F1Evaluation.F1Evaluation {
        input:
            tool1_label = tool1_label,
            tool2_label = tool2_label,
            additional_label = additional_label,
            signed_difference = signed_difference,
            roc_tables = roc_tables_3,
            preemptible = preemptible
    }

    # For the ROC plots it's ok to just flatten the array resulting from the scatter over i.
    Array[File] roc_plots_3 = flatten(roc_plots_2)

    # Also combine all plots into one image
    call MergePNGs as MergeFE {
        input:
            pngs = FEEvaluation.fe_plots,
            preemptible = preemptible
    }

    call MergePNGs as MergeF1 {
        input:
            pngs = F1Evaluation.f1_plots,
            preemptible = preemptible
    }

    call MergePNGs as MergeROC {
        input:
            pngs = roc_plots_3,
            preemptible = preemptible
    }

    output{
        Array[File] fe_plots = FEEvaluation.fe_plots
        Array[File] f1_plots = F1Evaluation.f1_plots
        Array[File] roc_plots = roc_plots_3
        File merged_fe_plots = MergeFE.plots
        File merged_f1_plots = MergeF1.plots
        File merged_roc_plots = MergeROC.plots
        File fe_summary = FEEvaluation.fe_summary
        File f1_summary = F1Evaluation.f1_summary
    }
}

task MergePNGs {
    input {
        Array[File] pngs
        Int? preemptible
    }

    command {
        convert -append ~{sep=" " pngs} plots.png
    }

    runtime {
        docker: "michaelgatzen/mergepngs"
        preemptible: select_first([preemptible, 0])
        memory: "8 GB"
        disks: "local-disk 20 HDD"
    }

    output{
        File plots = "plots.png"
    }
}

task RenameSummary {
    input {
        File input_summary
        String name
        Int? preemptible
    }

    command {
        mv ~{input_summary} ./~{name}.csv
    }

    runtime {
        docker: "ubuntu:20.04"
        preemptible: select_first([preemptible, 0])
        memory: "2 GB"
        disks: "local-disk 20 HDD"
    }

    output{
        File summary = "~{name}.csv"
    }
}