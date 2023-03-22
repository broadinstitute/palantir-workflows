version 1.0

import "LongReadRNABenchmarkTasks.wdl" as LongReadRNABenchmarkTasks

workflow LongReadRNABenchmark {
    input {
        File inputBAM
        File referenceGenome
        File reducedAnnotation
        String datasetName
        String dataType
        Int numThreads
    }

    call LongReadRNABenchmarkTasks.IsoQuant as IsoQuant {
        input:
            inputBAM = inputBAM,
            referenceGenome = referenceGenome,
            reducedAnnotation = reducedAnnotation,
            dataType = dataType,
            numThreads = numThreads
    }

    output {
        File isoQuantGTF = IsoQuant.isoQuantGTF
    }
}