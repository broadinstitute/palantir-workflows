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
#    call MethylationBenchmarkPerReferenceTasks.GenerateBWAIndex as GenerateBWAIndex {
#        input:
#            ref = ref
#    }

#    call MethylationBenchmarkReferenceTasks.GenerateFASTAIndex as GenerateFASTAIndex {
#        input:
#            ref = ref
#    }

#    call MethylationBenchmarkReferenceTasks.CreateSequenceDictionary as CreateSequenceDictionary {
#        input:
#            ref = ref
#    }

    output {
#        File bwaIdxAmb = GenerateBWAIndex.amb
#        File bwaIdxAnn = GenerateBWAIndex.ann
#        File bwaIdxBwt = GenerateBWAIndex.bwt
#        File bwaIdxPac = GenerateBWAIndex.pac
#        File bwaIdxSa = GenerateBWAIndex.sa
#        File fai = GenerateFASTAIndex.fai
#        File dict = CreateSequenceDictionary.dict
    }
}