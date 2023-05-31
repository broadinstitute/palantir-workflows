version 1.0

import "MethylationBenchmarkPerSampleTasks.wdl" as MethylationBenchmarkPerSampleTasks

workflow MethylationBenchmarkPerSample {
    input {
        String sampleId
        File ref
        File amb
        File ann
        File bwt
        File pac
        File sa
        File fq1gz
        File fq2gz
        File targets
    }

    call MethylationBenchmarkPerSampleTasks.DownsampleReads as DownsampleReadsFq1 {
        input:
            fqgz = fq1gz,
    }

    call MethylationBenchmarkPerSampleTasks.DownsampleReads as DownsampleReadsFq2 {
        input:
            fqgz = fq2gz,
    }

    call MethylationBenchmarkPerSampleTasks.TrimAdapters as TrimAdapters {
        input:
            fq1 = DownsampleReadsFq1.fqDownsampled,
            fq2 = DownsampleReadsFq2.fqDownsampled
    }

#    call MethylationBenchmarkPerSampleTasks.FastQC as FastQC {
#        input:
#            fq1gz = TrimAdapters.fq1Trimmed,
#            fq2gz = TrimAdapters.fq2Trimmed
#    }
#
#    call MethylationBenchmarkPerSampleTasks.BWAMethAlign as BWAMethAlign {
#        input:
#            sampleId = sampleId,
#            ref = ref,
#            amb = amb,
#            ann = ann,
#            bwt = bwt,
#            pac = pac,
#            sa = sa,
#            fq1 = TrimAdapters.fq1Trimmed,
#            fq2 = TrimAdapters.fq2Trimmed
#    }
#
#    call MethylationBenchmarkPerSampleTasks.SAMBamba as SAMBamba {
#        input:
#            ref = ref,
#            bam = BWAMethAlign.bam
#    }

    output {
        File fq1Downsampled = DownsampleReadsFq1.fqDownsampled
        File fq2Downsampled = DownsampleReadsFq2.fqDownsampled
        File fq1Trimmed = TrimAdapters.fq1Trimmed
        File fq2Trimmed = TrimAdapters.fq2Trimmed
#        File fastqcReportFq1 = FastQC.htmlReportFq1
#        File fastqcReportFq2 = FastQC.htmlReportFq2
#        File bam = BWAMethAlign.bam
#        File sortedBam = SAMBamba.sortedBam
    }
}