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

#    call LongReadRNABenchmarkTasks.Tama as Tama {
#        input:
#            inputBAM = inputBAM,
#            inputBAMIndex = inputBAMIndex,
#            referenceGenome = referenceGenome,
#            referenceGenomeIndex = referenceGenomeIndex,
#            datasetName = datasetName
#    }

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

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFIsoQuant {
        input:
            inputGTF = IsoQuant.isoQuantGTF,
            toolName = "isoquant"
    }

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFStringTie {
        input:
            inputGTF = StringTie.stringTieGTF,
            toolName = "stringtie"
    }

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFBambu {
        input:
            inputGTF = Bambu.bambuGTF,
            inputCounts = Bambu.bambuGTFCounts,
            toolName = "bambu"
    }

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFFlair {
        input:
            inputGTF = Flair.flairGTF,
            toolName = "flair"
    }

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFTalon {
        input:
            inputGTF = Talon.talonGTF,
            toolName = "talon"
    }

    call LongReadRNABenchmarkTasks.SplitGTF as SplitGTFFlames {
        input:
            inputGTF = Flames.flamesGFF,
            toolName = "flames"
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisIsoQuant {
        input:
            inputFullGTF = SplitGTFIsoQuant.full,
            inputKnownGTF = SplitGTFIsoQuant.known,
            inputNovelGTF = SplitGTFIsoQuant.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisStringTie {
        input:
            inputFullGTF = SplitGTFStringTie.full,
            inputKnownGTF = SplitGTFStringTie.known,
            inputNovelGTF = SplitGTFStringTie.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisBambu {
        input:
            inputFullGTF = SplitGTFBambu.full,
            inputKnownGTF = SplitGTFBambu.known,
            inputNovelGTF = SplitGTFBambu.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisFlair {
        input:
            inputFullGTF = SplitGTFFlair.full,
            inputKnownGTF = SplitGTFFlair.known,
            inputNovelGTF = SplitGTFFlair.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisTalon {
        input:
            inputFullGTF = SplitGTFTalon.full,
            inputKnownGTF = SplitGTFTalon.known,
            inputNovelGTF = SplitGTFTalon.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationAnalysis as ReducedAnnotationAnalysisFlames {
        input:
            inputFullGTF = SplitGTFFlames.full,
            inputKnownGTF = SplitGTFFlames.known,
            inputNovelGTF = SplitGTFFlames.novel,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis as ReferenceFreeAnalysisIsoQuant {
        input:
            inputGTF = IsoQuantReferenceFree.isoQuantGTF,
            expressedGTF = expressedGTF
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis as ReferenceFreeAnalysisStringTie {
        input:
            inputGTF = StringTieReferenceFree.stringTieGTF,
            expressedGTF = expressedGTF
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis as ReferenceFreeAnalysisIsoSeq {
        input:
            inputGTF = IsoSeq.isoSeqGFF,
            expressedGTF = expressedGTF
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis as ReferenceFreeAnalysisCupcake {
        input:
            inputGTF = Cupcake.cupcakeGFF,
            expressedGTF = expressedGTF
    }

#    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysis as ReferenceFreeAnalysisTama {
#        input:
#            inputGTF = Tama.tamaGTF,
#            expressedGTF = expressedGTF
#    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisIsoQuant {
        input:
            toolName = "isoquant",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisStringTie {
        input:
            toolName = "stringtie",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisBambu {
        input:
            toolName = "bambu",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisFlair {
        input:
            toolName = "flair",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisTalon {
        input:
            toolName = "talon",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoAnalysis as DenovoAnalysisFlames {
        input:
            toolName = "flames",
            gtfList = [IsoQuant.isoQuantGTF, StringTie.stringTieGTF, Bambu.bambuGTF, Flair.flairGTF, Talon.talonGTF, Flames.flamesGFF]
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsIsoQuant {
        input:
            trackingFile = DenovoAnalysisIsoQuant.tracking,
            toolName = "isoquant",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsStringTie {
        input:
            trackingFile = DenovoAnalysisStringTie.tracking,
            toolName = "stringtie",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsBambu {
        input:
            trackingFile = DenovoAnalysisBambu.tracking,
            toolName = "bambu",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsFlair {
        input:
            trackingFile = DenovoAnalysisFlair.tracking,
            toolName = "flair",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsTalon {
        input:
            trackingFile = DenovoAnalysisTalon.tracking,
            toolName = "talon",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.DenovoStats as DenovoStatsFlames {
        input:
            trackingFile = DenovoAnalysisFlames.tracking,
            toolName = "flames",
            numTools = 6
    }

    call LongReadRNABenchmarkTasks.SummarizeAnalysise as SummarizeAnalysisReduced {
        input:
            inputList = [ReducedAnnotationAnalysisIsoQuant.novel, ReducedAnnotationAnalysisStringTie.novel, ReducedAnnotationAnalysisBambu.novel ,ReducedAnnotationAnalysisFlair.novel ,ReducedAnnotationAnalysisTalon.novel, ReducedAnnotationAnalysisFlames.novel],
            toolNames = ["isoquant" , "stringtie", "bambu", "flair", "talon", "flames"],
            datasetName = datasetName,
            analysisType = "reduced"
    }

    call LongReadRNABenchmarkTasks.SummarizeAnalysise as SummarizeAnalysisReffree {
        input:
            inputList = [ReferenceFreeAnalysisIsoQuant.stats, ReferenceFreeAnalysisStringTie.stats, ReferenceFreeAnalysisIsoSeq.stats ,ReferenceFreeAnalysisCupcake.stats],
            toolNames = ["isoquant" , "stringtie", "isoseq", "cupcake"],
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
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantOut = IsoQuant.isoQuantOut
        File isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantGTF
        File isoQuantDenovoOut = IsoQuantReferenceFree.isoQuantOut
        File stringTieGTF = StringTie.stringTieGTF
        File stringTieDenovoGTF = StringTieReferenceFree.stringTieGTF
        File bambuGTF = Bambu.bambuGTF
        File bambuGTFCounts = Bambu.bambuGTFCounts
        File bambuOut = Bambu.bambuOut
        File flairGTF = Flair.flairGTF
        File talonGTF = Talon.talonGTF
        File isoSeqGFF = IsoSeq.isoSeqGFF
        #File tamaGTF = Tama.tamaGTF
        File flamesGFF = Flames.flamesGFF
        File cupcakeGFF = Cupcake.cupcakeGFF
        File isoQuantFull = SplitGTFIsoQuant.full
        File isoQuantKnown = SplitGTFIsoQuant.known
        File isoQuantNovel = SplitGTFIsoQuant.novel
        File isoQuantTracking = DenovoAnalysisIsoQuant.tracking
        File isoQuantReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisIsoQuant.novel
        File isoQuantReferenceFreeAnalysisStats = ReferenceFreeAnalysisIsoQuant.stats
        File stringTieFull = SplitGTFStringTie.full
        File stringTieKnown = SplitGTFStringTie.known
        File stringTieNovel = SplitGTFStringTie.novel
        File stringTieTracking = DenovoAnalysisStringTie.tracking
        File stringTieReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisStringTie.novel
        File stringTieReferenceFreeAnalysisStats = ReferenceFreeAnalysisStringTie.stats
        File bambuFull = SplitGTFBambu.full
        File bambuKnown = SplitGTFBambu.known
        File bambuNovel = SplitGTFBambu.novel
        File bambuTracking = DenovoAnalysisBambu.tracking
        File bambuReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisBambu.novel
        File flairFull = SplitGTFFlair.full
        File flairKnown = SplitGTFFlair.known
        File flairNovel = SplitGTFFlair.novel
        File flairTracking = DenovoAnalysisFlair.tracking
        File flairReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisFlair.novel
        File talonFull = SplitGTFTalon.full
        File talonKnown = SplitGTFTalon.known
        File talonNovel = SplitGTFTalon.novel
        File talonTracking = DenovoAnalysisTalon.tracking
        File talonReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisTalon.novel
        File flamesFull = SplitGTFFlames.full
        File flamesKnown = SplitGTFFlames.known
        File flamestNovel = SplitGTFFlames.novel
        File flamesTracking = DenovoAnalysisFlames.tracking
        File flamesReducedAnnotationAnalysisStats = ReducedAnnotationAnalysisFlames.novel
        File isoSeqReferenceFreeAnalysisStats = ReferenceFreeAnalysisIsoSeq.stats
        File cupcakeReferenceFreeAnalysisStats = ReferenceFreeAnalysisCupcake.stats
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        File isoQuantReferenceFreeMonitoringLog = IsoQuantReferenceFree.monitoringLog
        File stringTieMonitoringLog = StringTie.monitoringLog
        File stringTieReferenceFreeMonitoringLog = StringTieReferenceFree.monitoringLog
        File bambuMonitoringLog = Bambu.monitoringLog
        File flairMonitoringLog = Flair.monitoringLog
        File talonMonitoringLog = Talon.monitoringLog
        File isoSeqMonitoringLog = IsoSeq.monitoringLog
        #File tamaMonitoringLog = Tama.monitoringLog
        File flamesMonitoringLog = Flames.monitoringLog
        File cupcakeMonitoringLog = Cupcake.monitoringLog
        File reducedAnalysisSummary = SummarizeAnalysisReduced.summary
        File referenceFreeAnalysisSummary = SummarizeAnalysisReffree.summary
        File reducedAnalysisSummaryPlots = PlotAnalysisSummaryReduced.analysisSummaryPlots
        File referenceFreeAnalysisSummaryPlots = PlotAnalysisSummaryReffree.analysisSummaryPlots
    }
}