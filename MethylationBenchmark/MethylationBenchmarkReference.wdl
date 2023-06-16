version 1.0

import "MethylationBenchmarkReferenceTasks.wdl" as MethylationBenchmarkReferenceTasks

workflow MethylationBenchmarkReference {
    input {
        File ref
    }

    call MethylationBenchmarkReferenceTasks.GenerateBWAMethIndex{
        input:
            ref = ref
    }

    call MethylationBenchmarkReferenceTasks.GenerateFASTAIndex {
        input:
            ref = ref
    }

    call MethylationBenchmarkReferenceTasks.CreateSequenceDictionary {
        input:
            ref = ref
    }

    output {
        File bwamethIndex = GenerateBWAMethIndex.bwamethIndex
        File fastaIndex = GenerateFASTAIndex.fai
        File sequenceDictionary = CreateSequenceDictionary.dict
    }
}