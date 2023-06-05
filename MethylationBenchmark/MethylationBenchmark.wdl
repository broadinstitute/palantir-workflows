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

    call MethylationBenchmarkTasks.TrimAdapters as TrimAdapters {
        input:
            fq1 = fq1,
            fq2 = fq2
    }

    call MethylationBenchmarkTasks.FastQC as FastQC {
        input:
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed
    }

#    call MethylationBenchmarkTasks.BWAMethIndex as BWAMethIndex {
#        input:
#            ref = ref
#    }

    call MethylationBenchmarkTasks.BWAMethAlign as BWAMethAlign {
        input:
            sampleId = sampleId,
            ref = ref,
            fq1 = TrimAdapters.fq1Trimmed,
            fq2 = TrimAdapters.fq2Trimmed
    }

    call MethylationBenchmark.SAMBamba as SAMBamba {
        input:
            ref = ref,
            bam = BWAMethAlign.bam
    }

    output {
        #File fq1Downsampled = DownsampleReadsFq1.fqDownsampled
        #File fq2Downsampled = DownsampleReadsFq2.fqDownsampled
        File fq1Trimmed = TrimAdapters.fq1Trimmed
        File fq2Trimmed = TrimAdapters.fq2Trimmed
        File qcFq1 = FastQC.qcFq1
        File qcFq2 = FastQC.qcFq2
        File bam = BWAMethAlign.bam
        File sortedBam = SAMBamba.sortedBam
    }
}