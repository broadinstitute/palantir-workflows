version 1.0

import "LongReadRNABenchmarkTasks.wdl" as LongReadRNABenchmarkTasks

workflow LongReadRNABenchmark {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        File expressedGTF
        File expressedKeptGTF
        File excludedGTF
        String datasetName
        String dataType
    }

    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuant {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType
    }

    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuantReferenceFree {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName,
            dataType = dataType
    }

    call LongReadRNABenchmarkTasks.StringTie as StringTie {
        input:
            inputBAM = inputBAM,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.StringTie as StringTieReferenceFree {
        input:
            inputBAM = inputBAM,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.Bambu as Bambu {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType
    }

    call LongReadRNABenchmarkTasks.Flair as Flair {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.Talon as Talon {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType
    }

    call LongReadRNABenchmarkTasks.IsoSeq as IsoSeq {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.Flames as Flames {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.Cupcake as Cupcake {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            datasetName = datasetName
    }

    Array[File] gtfListReduced = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    Array[File] gtfListReferenceFree = [IsoQuantReferenceFree.isoQuantGTF, StringTieReferenceFree.stringTieGTF, IsoSeq.isoSeqGFF, Cupcake.cupcakeGFF]
    Array[String] toolNamesReduced = ["isoquant", "stringtie", "bambu", "flair", "talon", "flames"]
    Array[String] toolNamesReferenceFree = ["isoquant", "stringtie", "isoseq", "cupcake"]

    scatter(gtfAndTool in zip(gtfListReduced, toolNamesReduced)) {
        File gtf = gtfAndTool.left
        String tool = gtfAndTool.right

        if (tool != "bambu") {
            call LongReadRNABenchmarkTasks.SplitGTF {
                input:
                    inputGTF = gtf,
                    toolName = tool
            }
        }

        if (tool == "bambu") {
            call LongReadRNABenchmarkTasks.SplitGTF {
                input:
                    inputGTF = gtf,
                    inputCounts = Bambu.bambuCounts,
                    toolName = tool
            }
        }

        call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis {
            input:
                inputFullGTF = SplitGTF.full,
                inputKnownGTF = SplitGTF.known,
                inputNovelGTF = SplitGTF.novel,
                expressedGTF = expressedGTF,
                expressedKeptGTF = expressedKeptGTF,
                excludedGTF = excludedGTF
        }
    }

    scatter(gtf in gtfListReferenceFree) {
        call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis {
            input:
                inputGTF = gtf,
                expressedGTF = expressedGTF
        }
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisFull {
        input:
            splitType = "full",
            gtfList = SplitGTF.full
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisKnown {
        input:
            splitType = "known",
            gtfList = SplitGTF.known
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisNovel {
        input:
            splitType = "novel",
            gtfList = SplitGTF.novel
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsFull {
        input:
            splitType = "full",
            trackingFile = DenovoAnalysisFull.tracking,
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsKnown {
        input:
            splitType = "known",
            trackingFile = DenovoAnalysisKnown.tracking,
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsNovel {
        input:
            splitType = "novel",
            trackingFile = DenovoAnalysisNovel.tracking,
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.SummarizeAnalysis as SummarizeAnalysisReduced {
        input:
            inputList = ReducedAnnotationAnalysis.novel,
            toolNames = toolNamesReduced,
            datasetName = datasetName,
            analysisType = "reduced"
    }

    call LongReadRNABenchmarkTasks.SummarizeAnalysis as SummarizeAnalysisReffree {
        input:
            inputList = ReferenceFreeAnalysis.stats,
            toolNames = toolNamesReferenceFree,
            datasetName = datasetName,
            analysisType = "reffree"
    }

    call LongReadRNABenchmarkTasks.PlotAnalysisSummary as PlotAnalysisSummaryReduced {
        input:
            summary = SummarizeAnalysisReduced.summary,
            datasetName = datasetName,
            analysisType = "reduced"
    }

    call LongReadRNABenchmarkTasks.PlotAnalysisSummary as PlotAnalysisSummaryReffree {
        input:
            summary = SummarizeAnalysisReffree.summary,
            datasetName = datasetName,
            analysisType = "reffree"
    }

    output {
        File reducedAnalysisSummary = SummarizeAnalysisReduced.summary
        File referenceFreeAnalysisSummary = SummarizeAnalysisReffree.summary
        File reducedAnalysisSummaryPlots = PlotAnalysisSummaryReduced.analysisSummaryPlots
        File referenceFreeAnalysisSummaryPlots = PlotAnalysisSummaryReffree.analysisSummaryPlots
    }
}