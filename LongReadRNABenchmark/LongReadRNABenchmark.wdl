version 1.0

import "IsoQuant.wdl" as IsoQuantWorkflow
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

    call IsoQuantWorkflow.IsoQuant as IsoQuant {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType
    }

    call IsoQuantWorkflow.IsoQuant as IsoQuantReferenceFree {
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

    # Note: Make sure that your toolNames arrays match the order of your gtfList arrays.
    # If they don't match, you may not get an error but you will get incorrect results.
    Array[File] gtfListReduced = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    Array[File] gtfListReferenceFree = [IsoQuantReferenceFree.isoQuantGTF, StringTieReferenceFree.stringTieGTF, IsoSeq.isoSeqGFF, Cupcake.cupcakeGFF]
    Array[String] toolNamesReduced = ["isoquant", "stringtie", "bambu", "flair", "talon", "flames"]
    Array[String] toolNamesReferenceFree = ["isoquant", "stringtie", "isoseq", "cupcake"]

    scatter(gtfAndTool in zip(gtfListReduced, toolNamesReduced)) {
        File gtf = gtfAndTool.left
        String tool = gtfAndTool.right

        call LongReadRNABenchmarkTasks.SplitGTF {
            input:
                inputGTF = gtf,
                inputCounts = Bambu.bambuCounts,
                toolName = tool
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
            datasetName = datasetName,
            splitType = "full",
            gtfList = SplitGTF.full
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisKnown {
        input:
            datasetName = datasetName,
            splitType = "known",
            gtfList = SplitGTF.known
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisNovel {
        input:
            datasetName = datasetName,
            splitType = "novel",
            gtfList = SplitGTF.novel
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsFull {
        input:
            datasetName = datasetName,
            splitType = "full",
            trackingFile = DenovoAnalysisFull.tracking,
            numTools = 6,
            toolNames = toolNamesReduced
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsKnown {
        input:
            datasetName = datasetName,
            splitType = "known",
            trackingFile = DenovoAnalysisKnown.tracking,
            numTools = 6,
            toolNames = toolNamesReduced
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsNovel {
        input:
            datasetName = datasetName,
            splitType = "novel",
            trackingFile = DenovoAnalysisNovel.tracking,
            numTools = 6,
            toolNames = toolNamesReduced
    }

    call LongReadRNABenchmarkTasks.PlotDenovoStats as PlotDenovoStatsFull {
        input:
            stats = DenovoStatsFull.denovoStats,
            datasetName = datasetName,
            splitType = "full"
    }

    call LongReadRNABenchmarkTasks.PlotDenovoStats as PlotDenovoStatsKnown {
        input:
            stats = DenovoStatsKnown.denovoStats,
            datasetName = datasetName,
            splitType = "known"
    }

    call LongReadRNABenchmarkTasks.PlotDenovoStats as PlotDenovoStatsNovel {
        input:
            stats = DenovoStatsNovel.denovoStats,
            datasetName = datasetName,
            splitType = "novel"
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

    call LongReadRNABenchmarkTasks.GenerateSplitFreeTracking as GenerateSplitFreeTrackingIsoQuant {
        input:
            datasetName = datasetName,
            toolGTF = IsoQuant.isoQuantGTF,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF
    }

    call LongReadRNABenchmarkTasks.SplitFreeStats as SplitFreeStatsIsoQuant {
        input:
            trackingFile = GenerateSplitFreeTrackingIsoQuant.tracking,
            toolName = "isoquant",
            datasetName = datasetName
    }

    output {
        File reducedAnalysisSummary = SummarizeAnalysisReduced.summary
        File referenceFreeAnalysisSummary = SummarizeAnalysisReffree.summary
        File reducedAnalysisSummaryPlots = PlotAnalysisSummaryReduced.analysisSummaryPlots
        File referenceFreeAnalysisSummaryPlots = PlotAnalysisSummaryReffree.analysisSummaryPlots
        File fullDenovoStatsPlot = PlotDenovoStatsFull.denovoStatsPlot
        File knownDenovoStatsPlot = PlotDenovoStatsKnown.denovoStatsPlot
        File novelDenovoStatsPlot = PlotDenovoStatsNovel.denovoStatsPlot
    }
}