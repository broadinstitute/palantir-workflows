version 1.0

import "MethylationBenchmarkPerReferenceTasks.wdl" as MethylationBenchmarkPerReferenceTasks

workflow MethylationBenchmark {
    input {
        File ref
    }

    call MethylationBenchmarkPerReferenceTasks.GenerateBWAIndex as GenerateBWAIndex {
        input:
            ref = ref
    }

    call MethylationBenchmarkPerReferenceTasks.GenerateFASTAIndex as GenerateFASTAIndex {
        input:
            ref = ref
    }

    call MethylationBenchmarkPerReferenceTasks.CreateSequenceDictionary as CreateSequenceDictionary {
        input:
            ref = ref
    }

    output {
        File bwaIdxAmb = GenerateBWAIndex.amb
        File bwaIdxAnn = GenerateBWAIndex.ann
        File bwaIdxBwt = GenerateBWAIndex.bwt
        File bwaIdxPac = GenerateBWAIndex.pac
        File bwaIdxSa = GenerateBWAIndex.sa
        File fai = GenerateFASTAIndex.fai
        File dict = CreateSequenceDictionary.dict
    }
}