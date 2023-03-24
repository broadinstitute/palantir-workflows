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

    #call LongReadRNABenchmarkTasks.StringTie as StringTie {
    #    input:
    #        inputBAM = inputBAM,
    #        referenceAnnotation = referenceAnnotation,
    #        datasetName = datasetName,
    #        numThreads = numThreads
    #}

    output {
        File isoQuantGTF = IsoQuant.isoQuantGTF
        File isoQuantDB = IsoQuant.isoQuantDB
        File isoQuantMonitoringLog = IsoQuant.monitoringLog
        #File stringTieGTF = StringTie.stringTieGTF
    }
}