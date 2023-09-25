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

    call EMSeqTasks.Fastp {
        input:
            sampleId = sampleId,
            fq1 = fq1,
            fq2 = fq2
    }

    call EMSeqTasks.Bwameth {
        input:
            sampleId = sampleId,
            fq1 = Fastp.filteredFq1,
            fq2 = Fastp.filteredFq2,
            ref = ref,
            refIdx = refIdx,
            bwamethIdx = bwamethIdx
    }

    call EMSeqTasks.MarkNonconvertedReads {
        input:
            sampleId = sampleId,
            sam = Bwameth.sam,
            ref = ref,
            refIdx = refIdx,
    }

    call EMSeqTasks.Sambamba {
        input:
            sampleId = sampleId,
            sam = MarkNonconvertedReads.ncMarkedSam
    }

    call EMSeqTasks.MarkDuplicates {
        input:
            sampleId = sampleId,
            bam = Sambamba.bam
    }

    call EMSeqTasks.MethylDackelMbias {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.mdBam,
            bai = MarkDuplicates.mdBai,
            ref = ref,
            refIdx = refIdx
    }

    call EMSeqTasks.MethylDackelExtract {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.mdBam,
            bai = MarkDuplicates.mdBai,
            ref = ref,
            refIdx = refIdx
    }

    output {
        File filteredFq1 = Fastp.filteredFq1
        File filteredFq2 = Fastp.filteredFq2
        File fastpReport = Fastp.fastpReport
        File sam = Bwameth.sam
        File ncMarkedSam = MarkNonconvertedReads.ncMarkedSam
        File nonconvertedReadCounts = MarkNonconvertedReads.nonconvertedReadCounts
        File sambambaBam = Sambamba.bam
        File bam = MarkDuplicates.mdBam
        File bai = MarkDuplicates.mdBai
        File samblasterLog = MarkDuplicates.samblasterLog
        File combinedMbias = MethylDackelMbias.combinedMbias
        File chnSvgOB = MethylDackelMbias.chnSvgOB
        File chnSvgOT = MethylDackelMbias.chnSvgOT
        File cpgSvgOB = MethylDackelMbias.cpgSvgOB
        File cpgSvgOT = MethylDackelMbias.cpgSvgOT
        #File methylKitCHG = MethylDackelExtract.methylKitCHG
        #File methylKitCHH = MethylDackelExtract.methylKitCHH
        #File methylKitCpG = MethylDackelExtract.methylKitCpG
        File cpgBedGraph = MethylDackelExtract.cpgBedGraph
    }
}