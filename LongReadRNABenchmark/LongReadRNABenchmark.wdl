version 1.0

import "LongReadRNABenchmarkTasks.wdl" as LongReadRNABenchmarkTasks

workflow LongReadRNABenchmark {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
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

    call LongReadRNABenchmarkTasks.ReducedAnnotationGFFCompare as ReducedAnnotationGFFCompare {
        input:
            reducedAnnotationDB = IsoQuant.isoQuantDB,
            isoQuantGTF = IsoQuant.isoQuantGTF,
            stringTieGTF = StringTie.stringTieGTF,
            datasetName = datasetName
    }

    output {
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantDenovoGTF = IsoQuantReferenceFree.isoQuantDenovoGTF
        File stringTieGTF = StringTie.stringTieGTF
        File stringTieDenovoGTF = StringTieReferenceFree.stringTieDenovoGTF
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        File isoQuantReferenceFreeMonitoringLog = IsoQuantReferenceFree.monitoringLog
        File stringTieMonitoringLog = StringTie.monitoringLog
        File stringTieReferenceFreeMonitoringLog = StringTieReferenceFree.monitoringLog
    }
}