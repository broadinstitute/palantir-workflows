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
        Int numThreads
    }

    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuant {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.IsoQuantReferenceFree as IsoQuantReferenceFree {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName,
            dataType = dataType,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.StringTie as StringTie {
        input:
            inputBAM = inputBAM,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.StringTieReferenceFree as StringTieReferenceFree {
        input:
            inputBAM = inputBAM,
            datasetName = datasetName,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.Bambu as Bambu {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.Flair as Flair {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.Talon as Talon {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType,
            numThreads = numThreads
    }

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompare {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            expressedGTF = expressedGTF,
            expressedKeptGTF = expressedKeptGTF,
            excludedGTF = excludedGTF,
            isoQuantGTF = IsoQuant.isoQuantGTF,
            stringTieGTF = StringTie.stringTieGTF,
            bambuGTF = Bambu.bambuGTF,
            bambuGTFCounts = Bambu.bambuGTFCounts,
            flairGTF = Flair.flairGTF,
            talonGTF = Talon.talonGTF,
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
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReferenceFreeGFFCompare as ReferenceFreeGFFCompare {
        input:
            isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF,
            stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF,
            expressedGTF = expressedGTF,
            datasetName = datasetName
    }

    call LongReadRNABenchmarkTasks.ReducedAnalysisSummarize as ReducedAnalysisSummarize {
        input:
            reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompare.gffCompareOutputIsoQuant,
            reducedGffCompareOutStringTie = ReducedAnnotationGFFCompare.gffCompareOutputStringTie,
            reducedGffCompareOutBambu = ReducedAnnotationGFFCompare.gffCompareOutputBambu,
            reducedGffCompareOutFlair = ReducedAnnotationGFFCompare.gffCompareOutputFlair,
            reducedGffCompareOutTalon = ReducedAnnotationGFFCompare.gffCompareOutputTalon,
            datasetName = datasetName
    }

    output {
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantOut = IsoQuant.isoQuantOut
        File isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF
        File isoQuantDenovoOut = IsoQuantReferenceFree.isoQuantDenovoOut
        File stringTieGTF = StringTie.stringTieGTF
        File stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF
        File bambuGTF = Bambu.bambuGTF
        File bambuGTFCounts = Bambu.bambuGTFCounts
        File bambuOut = Bambu.bambuOut
        File flairGTF = Flair.flairGTF
        File talonGTF = Talon.talonGTF
        File denovoAnnotationGFFCompareOut = DenovoAnnotationGFFCompare.gffCompareOutput
        File reducedGffCompareOutIsoQuant = ReducedAnnotationGFFCompare.gffCompareOutputIsoQuant
        File reducedGffCompareOutStringTie = ReducedAnnotationGFFCompare.gffCompareOutputStringTie
        File reducedGffCompareOutBambu = ReducedAnnotationGFFCompare.gffCompareOutputBambu
        File reducedGffCompareOutFlair = ReducedAnnotationGFFCompare.gffCompareOutputFlair
        File reducedGffCompareOutTalon = ReducedAnnotationGFFCompare.gffCompareOutputTalon
        File referenceFreeGFFCompareOut = ReferenceFreeGFFCompare.gffCompareOutput
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        File isoQuantReferenceFreeMonitoringLog = IsoQuantReferenceFree.monitoringLog
        File stringTieMonitoringLog = StringTie.monitoringLog
        File stringTieReferenceFreeMonitoringLog = StringTieReferenceFree.monitoringLog
        File bambuMonitoringLog = Bambu.monitoringLog
        File flairMonitoringLog = Flair.monitoringLog
        File talonMonitoringLog = Talon.monitoringLog
        File reducedAnnotationGFFCompareMonitoringLog = ReducedAnnotationGFFCompare.monitoringLog
        File denovoAnnotationGFFCompareMonitoringLog = DenovoAnnotationGFFCompare.monitoringLog
        File referenceFreeGFFCompareMonitoringLog = ReferenceFreeGFFCompare.monitoringLog
        File reducedAnalysisSummary = ReducedAnalysisSummarize.reducedAnalysisSummary
        File reducedAnalysisAccuracyPlots = ReducedAnalysisSummarize.reducedAnalysisAccuracyPlots
    }
}