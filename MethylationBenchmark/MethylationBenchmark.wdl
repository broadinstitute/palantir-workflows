version 1.0

import "MethylationBenchmarkTasks.wdl" as MethylationBenchmarkTasks

workflow MethylationBenchmark {
    input {
        String sampleId
        File fq1
        File fq2
        File ref
        File refIdx
        File targets
        #File amb
        #File ann
        #File bwt
        #File pac
        #File sa
    }

#    call MethylationBenchmarkPerSampleTasks.DownsampleReads as DownsampleReadsFq1 {
#        input:
#            fqgz = fq1gz,
#    }
#
#    call MethylationBenchmarkPerSampleTasks.DownsampleReads as DownsampleReadsFq2 {
#        input:
#            fqgz = fq2gz,
#    }

    call MethylationBenchmarkTasks.TrimAdapters {
        input:
            fq1 = fq1,
            fq2 = fq2
    }

    call MethylationBenchmarkTasks.FastQC {
        input:
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed
    }

#    call MethylationBenchmarkTasks.BWAMethIndex as BWAMethIndex {
#        input:
#            ref = ref
#    }

    call MethylationBenchmarkTasks.BWAMethAlign {
        input:
            sampleId = sampleId,
            ref = ref,
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed
    }

    call MethylationBenchmarkTasks.SAMBambaFilter {
        input:
            ref = ref,
            sam = BWAMethAlign.sam
    }

    call MethylationBenchmarkTasks.SAMBambaSort {
        input:
            ref = ref,
            bam = SAMBambaFilter.filteredBam
    }

    call MethylationBenchmarkTasks.SamtoolsIndex as SamtoolsIndexSortedBam {
        input:
            bam = SAMBambaSort.sortedBam
    }

    call MethylationBenchmarkTasks.MarkDuplicates {
        input:
            sampleId = sampleId,
            bam = SAMBambaSort.sortedBam,
            ref = ref
    }

    call MethylationBenchmarkTasks.SamtoolsIndex as SamtoolsIndexMarkdupBam {
        input:
            bam = MarkDuplicates.markdupBam
    }

    call MethylationBenchmarkTasks.CreateSequenceDictionary {
        input:
            ref = ref
    }

    call MethylationBenchmarkTasks.BedToIntervalList {
        input:
            bed = targets,
            dict = CreateSequenceDictionary.dict
    }

    call MethylationBenchmarkTasks.CollectHsMetrics {
        input:
            sampleId = sampleId,
            ref = ref,
            refIdx = refIdx,
            bam = MarkDuplicates.markdupBam,
            intervals = BedToIntervalList.intervalList
    }

    call MethylationBenchmarkTasks.CollectMultipleMetrics {
        input:
            sampleId = sampleId,
            ref = ref,
            bam = MarkDuplicates.markdupBam
    }

    call MethylationBenchmarkTasks.MethylDackelMbias {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.markdupBam,
            bai = SamtoolsIndexMarkdupBam.bai,
            ref = ref,
            refIdx = refIdx
    }

    call MethylationBenchmarkTasks.MethylDackelCallCpG {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.markdupBam,
            bai = SamtoolsIndexMarkdupBam.bai,
            ref = ref,
            refIdx = refIdx,
            mbiasParams = MethylDackelMbias.mbiasParams
    }

    call MethylationBenchmarkTasks.CreateMoreSignificantFiguresForPercentMethylation {
        input:
            cgpBedGraph = MethylDackelCallCpG.cpgBedGraph
    }

    call MethylationBenchmarkTasks.MethylDackelGenerateCytosineReport {
        input:
            sampleId = sampleId,
            bam = MarkDuplicates.markdupBam,
            bai = SamtoolsIndexMarkdupBam.bai,
            ref = ref,
            refIdx = refIdx,
            mbiasParams = MethylDackelMbias.mbiasParams
    }

    output {
        #File fq1Downsampled = DownsampleReadsFq1.fqDownsampled
        #File fq2Downsampled = DownsampleReadsFq2.fqDownsampled
        File fq1Trimmed = TrimAdapters.fq1Trimmed
        File fq2Trimmed = TrimAdapters.fq2Trimmed
        File qcFq1 = FastQC.qcFq1
        File qcFq2 = FastQC.qcFq2
        File sam = BWAMethAlign.sam
        File filteredBam = SAMBambaFilter.filteredBam
        File sortedBam = SAMBambaSort.sortedBam
        File sortedBai = SamtoolsIndexSortedBam.bai
        File markdupBam = MarkDuplicates.markdupBam
        File markdupBai = SamtoolsIndexMarkdupBam.bai
        File markdupMetrics = MarkDuplicates.markdupMetrics
        File dict = CreateSequenceDictionary.dict
        File intervalList = BedToIntervalList.intervalList
        File hsMetrics = CollectHsMetrics.hsMetrics
        File perTargetCoverage = CollectHsMetrics.perTargetCoverage
        File gcBiasDetail = CollectMultipleMetrics.gcBiasDetail
        File gcBiasSummary = CollectMultipleMetrics.gcBiasSummary
        File gcBiasPdf = CollectMultipleMetrics.gcBiasPdf
        File insertSizeMetrics = CollectMultipleMetrics.insertSizeMetrics
        File insertSizeHistogram = CollectMultipleMetrics.insertSizeHistogram
        File OB = MethylDackelMbias.OB
        File OT = MethylDackelMbias.OT
        File mbiasParams = MethylDackelMbias.mbiasParams
        File cpgBedGraph = MethylDackelCallCpG.cpgBedGraph
        File processedCpGBedGraph = CreateMoreSignificantFiguresForPercentMethylation.processedCpGBedGraph
        File cytosineReport = MethylDackelGenerateCytosineReport.cytosineReport
    }
}