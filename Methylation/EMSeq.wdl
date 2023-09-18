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

    call EMSeqTasks.MarkDuplicates {
        input:
            sampleId = sampleId,
            bam = Mapping.bam
    }

    output {
        File bam = Mapping.bam
        File fastpReport = Mapping.fastpReport
        File nonconvertedReadCounts = Mapping.nonconvertedReadCounts
    }
}