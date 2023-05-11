version 1.0

import "MethylationBenchmarkTasks.wdl" as MethylationBenchmarkTasks

workflow MethylationBenchmark {
    input {
        File ref
        File fq1gz
        File fq2gz
        File targets
    }

    call MethylationBenchmarkTasks.GenerateBWAIndex as GenerateBWAIndex {
        input:
            ref = ref
    }

    call MethylationBenchmarkTasks.GenerateFASTAIndex as GenerateFASTAIndex {
        input:
            ref = ref
    }

    call MethylationBenchmarkTasks.CreateSequenceDictionary as CreateSequenceDictionary {
        input:
            ref = ref
    }

    #call MethylationBenchmarkTasks.FastQC as FastQC {
    #    input:
    #        fq1 = fq1,
    #        fq2 = fq2
    #}

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