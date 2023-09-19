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

    call EMSeqTasks.MethylDackelMbias {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.mdBam,
            bai = MarkDuplicates.mdBai,
            ref = ref,
            refIdx = refIdx
    }

    output {
        File bam = MarkDuplicates.mdBam
        File bai = MarkDuplicates.mdBai
        File fastpReport = Mapping.fastpReport
        File nonconvertedReadCounts = Mapping.nonconvertedReadCounts
        File samblasterLog = MarkDuplicates.samblasterLog
        File combinedMbias = MethylDackelMbias.combinedMbias
        File chnSvgOB = MethylDackelMbias.chnSvgOB
        File chnSvgOT = MethylDackelMbias.chnSvgOT
        File cpgSvgOB = MethylDackelMbias.cpgSvgOB
        File cpgSvgOT = MethylDackelMbias.cpgSvgOT
    }
}