version 1.0

import "EMSeqTasks.wdl" as EMSeqTasks

workflow EMSeq {
    input {
        String sampleId
        File fq1
        File fq2
        File ref
        File refIdx
        File bwamethIdx
    }

    call EMSeqTasks.Mapping {
        input:
            sampleId = sampleId,
            fq1 = fq1,
            fq2 = fq2,
            ref = ref,
            refIdx = refIdx,
            bwamethIdx = bwamethIdx
    }

    output {
        #File sam = Mapping.sam
    }
}