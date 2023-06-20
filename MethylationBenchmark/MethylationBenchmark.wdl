version 1.0

import "MethylationBenchmarkTasks.wdl" as MethylationBenchmarkTasks

workflow MethylationBenchmark {
    input {
        String sampleId
        File fq1
        File fq2
        File ref
        File refIdx
        File bwamethIdx
        File targets
        Boolean downsample
    }

    if(downsample)
    {
        call MethylationBenchmarkTasks.DownsampleReads as DownsampleReadsFq1 {
            input:
                fq = fq1,
        }

        call MethylationBenchmarkTasks.DownsampleReads as DownsampleReadsFq2 {
            input:
                fq = fq2,
        }
    }

    call MethylationBenchmarkTasks.TrimAdapters {
        input:
            fq1 = select_first([DownsampleReadsFq1.fqDownsampled, fq1]),
            fq2 = select_first([DownsampleReadsFq2.fqDownsampled, fq2])
    }

    call MethylationBenchmarkTasks.FastQC {
        input:
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed
    }

    call MethylationBenchmarkTasks.BWAMethAlign {
        input:
            sampleId = sampleId,
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed,
            ref = ref,
            refIdx = refIdx,
            bwamethIdx = bwamethIdx
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
            ref = ref,
            refIdx = refIdx
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
            refIdx = refIdx,
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
            cpgBedGraph = MethylDackelCallCpG.cpgBedGraph
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

    call MethylationBenchmarkTasks.CollectMethylationStatistics {
        input:
            sampleId = sampleId,
            originalSam = BWAMethAlign.sam,
            filteredBam = SAMBambaSort.sortedBam,
            filteredBai = SamtoolsIndexSortedBam.bai,
            cytosineReport = MethylDackelGenerateCytosineReport.cytosineReport
    }

    output {
        File? fq1Downsampled = DownsampleReadsFq1.fqDownsampled
        File? fq2Downsampled = DownsampleReadsFq2.fqDownsampled
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
        Float mappingEfficiency = CollectMethylationStatistics.mappingEfficiency
    }
}