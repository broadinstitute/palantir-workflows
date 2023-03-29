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
            # Talon does not like reference genomes with the ">chr13 13" format.
            # For now, input a different, processed version unlike the other tools
            #referenceGenome = referenceGenome,
            #referenceGenomeIndex = referenceGenomeIndex,
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

    output {
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF
        File stringTieGTF = StringTie.stringTieGTF
        File stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF
        File bambuGTF = Bambu.bambuGTF
        File bambuGTFCounts = Bambu.bambuGTFCounts
        File flairGTF = Flair.flairGTF
        File talonGTF = Talon.talonGTF
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        File isoQuantReferenceFreeMonitoringLog = IsoQuantReferenceFree.monitoringLog
        File stringTieMonitoringLog = StringTie.monitoringLog
        File stringTieReferenceFreeMonitoringLog = StringTieReferenceFree.monitoringLog
        File bambuMonitoringLog = Bambu.monitoringLog
        File flairMonitoringLog = Flair.monitoringLog
        File talonMonitoringLog = Talon.monitoringLog
    }
}