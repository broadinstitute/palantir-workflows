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

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareIsoQuant {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = IsoQuant.isoQuantGTF,
            toolName = "isoquant",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareStringTie {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = StringTie.stringTieGTF,
            toolName = "stringtie",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareBambu {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            counts = Bambu.bambuGTFCounts,
            inputGTF = Bambu.bambuGTF,
            toolName = "bambu",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareFlair {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Flair.flairGTF,
            toolName = "flair",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareTalon {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Talon.talonGTF,
            toolName = "talon",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompareFlames {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            inputGTF = Flames.flamesGFF,
            toolName = "flames",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.DenovoAnnotationGFFCompare as DenovoAnnotationGFFCompare {
        input:
            isoQuantGTF = IsoQuant.isoQuantGTF,
            stringTieGTF = StringTie.stringTieGTF,
            bambuGTF = Bambu.bambuGTF,
            bambuGTFCounts = Bambu.bambuGTFCounts,
            flairGTF = Flair.flairGTF,
            talonGTF = Talon.talonGTF,
            flamesGFF = Flames.flamesGFF,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompareIsoQuant {
        input:
            inputGTF = IsoQuantReferenceFree.isoQuantGTF,
            expressedGTF = expressedGTF,
            toolName = "isoquant",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompareStringTie {
        input:
            inputGTF = StringTieReferenceFree.stringTieGTF,
            expressedGTF = expressedGTF,
            toolName = "stringtie",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompareIsoSeq {
        input:
            inputGTF = IsoSeq.isoSeqGFF,
            expressedGTF = expressedGTF,
            toolName = "isoseq",
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompareCupcake {
        input:
            inputGTF = Cupcake.cupcakeGFF,
            expressedGTF = expressedGTF,
            toolName = "cupcake",
            datasetName = datasetName
    }

#    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompareTama {
#        input:
#            inputGTF = Tama.tamaGTF,
#            expressedGTF = expressedGTF,
#            toolName = "tama",
#            datasetName = datasetName
#    }

    call LongReadRNABenchmarkTasks.ReducedAnalysisSummarize as ReducedAnalysisSummarize {
        input:
            reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompareIsoQuant.gffCompareOutput,
            reducedGffCompareOutStringTie = ReducedAnnotationGFFCompareStringTie.gffCompareOutput,
            reducedGffCompareOutBambu = ReducedAnnotationGFFCompareBambu.gffCompareOutput,
            reducedGffCompareOutFlair = ReducedAnnotationGFFCompareFlair.gffCompareOutput,
            reducedGffCompareOutTalon = ReducedAnnotationGFFCompareTalon.gffCompareOutput,
            reducedGffCompareOutFlames = ReducedAnnotationGFFCompareFlames.gffCompareOutput,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeAnalysisSummarize as ReferenceFreeAnalysisSummarize {
        input:
            referenceFreeGffCompareOutIsoQuant = ReferenceFreeGFFCompareIsoQuant.gffCompareOutput,
            referenceFreeGffCompareOutStringTie = ReferenceFreeGFFCompareStringTie.gffCompareOutput,
            referenceFreeGffCompareOutIsoSeq = ReferenceFreeGFFCompareIsoSeq.gffCompareOutput,
            #referenceFreeGffCompareOutTama = ReferenceFreeGFFCompareTama.gffCompareOutput,
            referenceFreeGffCompareOutCupcake = ReferenceFreeGFFCompareCupcake.gffCompareOutput,
            datasetName = datasetName
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
        File denovoAnnotationGFFCompareOut = DenovoAnnotationGFFCompare.gffCompareOutput
        File reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompareIsoQuant.gffCompareOutput
        File reducedGffCompareOutStringTie = ReducedAnnotationGFFCompareStringTie.gffCompareOutput
        File reducedGffCompareOutBambu = ReducedAnnotationGFFCompareBambu.gffCompareOutput
        File reducedGffCompareOutFlair = ReducedAnnotationGFFCompareFlair.gffCompareOutput
        File reducedGffCompareOutTalon = ReducedAnnotationGFFCompareTalon.gffCompareOutput
        File reducedGffCompareOutFlames = ReducedAnnotationGFFCompareFlames.gffCompareOutput
        File referenceFreeGFFCompareOutIsoQuant = ReferenceFreeGFFCompareIsoQuant.gffCompareOutput
        File referenceFreeGFFCompareOutStringTie = ReferenceFreeGFFCompareStringTie.gffCompareOutput
        File referenceFreeGFFCompareOutIsoSeq = ReferenceFreeGFFCompareIsoSeq.gffCompareOutput
        File referenceFreeGFFCompareOutCupcake = ReferenceFreeGFFCompareCupcake.gffCompareOutput
        #File referenceFreeGFFCompareOutTama = ReferenceFreeGFFCompareTama.gffCompareOutput
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
        File reducedAnalysisSummary = ReducedAnalysisSummarize.reducedAnalysisSummary
        File reducedAnalysisAccuracyPlots = ReducedAnalysisSummarize.reducedAnalysisAccuracyPlots
        File referenceFreeAnalysisSummary = ReferenceFreeAnalysisSummarize.referenceFreeAnalysisSummary
        File referenceFreeAnalysisAccuracyPlots = ReferenceFreeAnalysisSummarize.referenceFreeAnalysisAccuracyPlots
    }
}