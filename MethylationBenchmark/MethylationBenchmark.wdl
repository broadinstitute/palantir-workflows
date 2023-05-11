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

    #call MethylationBenchmarkTasks.SamtoolsFaidx as SamtoolsFaidx {
    #    input:
    #        reference = reference
    #}

    #call MethylationBenchmarkTasks.CreateSequenceDictionary as CreateSequenceDictionary {
    #    input:
    #        reference = reference
    #}

    #call MethylationBenchmarkTasks.FastQC as FastQC {
    #    input:
    #        fq1 = fq1,
    #        fq2 = fq2
    #}

    output {
        #File amb = BWAIndex.amb
        #File ann = BWAIndex.ann
        #File bwt = BWAIndex.bwt
        #File pac = BWAIndex.pac
        #File sa = BWAIndex.sa
        #File fai = SamtoolsFaidx.fai
        #File dict = CreateSequenceDictionary.dict
    }
}