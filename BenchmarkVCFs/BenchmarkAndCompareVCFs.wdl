version 1.0

import "BenchmarkVCFs.wdl" as BenchmarkVCFs
import "CompareBenchmarks.wdl" as CompareBenchmarks

workflow BenchmarkAndCompareVCFs {
    input {
        Array[String] sample_ids
        Array[String] configurations

        Boolean include_counts = false

        Array[String]? order_of_samples
        Array[String]? order_of_configurations
        Array[Int]? deltas

        Int? preemptible
        # Variables on a input sample level
        Array[File] evalVcf
        Array[File] evalVcfIndex
        Array[File] truthVcf
        Array[File] confidenceInterval
        Array[File] truthVcfIndex

        # Set on an analysis level
        String? analysisRegion
        Array[String]? jexlVariantSelectors
        Array[String]? variantSelectorLabels
        File? gatkJarForAnnotation
        Array[String]? annotationNames
        String? vcfScoreField
        String? dummyInputForTerraCallCaching
        File reference
        File refIndex
        File refDict
        File hapMap
        Array[File] stratIntervals = []
        Array[String] stratLabels = []
        String referenceVersion
        Int? threadsVcfEval=2
        Boolean doIndelLengthStratification=true
        Int? preemptible
        String gatkTag="4.0.11.0"
        Boolean requireMatchingGenotypes=true
        Boolean enableRefOverlap = false
        Boolean passingOnly=true
    }


    scatter(i in range(length(evalVcf))) {
        call BenchmarkVCFs.Benchmark as Benchmark {
            input:
                evalVcf = evalVcf[i],
                evalLabel = sample_ids[i],
                evalVcfIndex = evalVcfIndex[i],
                truthVcf = truthVcf[i],
                confidenceInterval = confidenceInterval[i],
                truthLabel = "truth",
                truthVcfIndex = truthVcfIndex[i],

                analysisRegion = analysisRegion,
                reference = reference,
                refIndex = refIndex,
                refDict = refDict,
                hapMap = hapMap,
                stratIntervals = stratIntervals,
                stratLabels = stratLabels,
                jexlVariantSelectors = jexlVariantSelectors,
                variantSelectorLabels = variantSelectorLabels,
                referenceVersion = referenceVersion,
                threadsVcfEval = threadsVcfEval,
                doIndelLengthStratification = doIndelLengthStratification,
                preemptible = preemptible,
                gatkTag = gatkTag,
                requireMatchingGenotypes = requireMatchingGenotypes,
                gatkJarForAnnotation = gatkJarForAnnotation,
                annotationNames = annotationNames,
                enableRefOverlap = enableRefOverlap,
                passingOnly = passingOnly,
                vcfScoreField = vcfScoreField,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }
    }

    call CompareBenchmarks.CompareBenchmarks as CompareBenchmarks {
        input:
            benchmark_summaries = Benchmark.summary,
            sample_ids = sample_ids,
            configurations = configurations,
            stratifiers = stratLabels,
            include_counts = include_counts,
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
        Array[File] gc_plots = CompareBenchmarks.gc_plots
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