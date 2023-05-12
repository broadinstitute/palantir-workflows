version 1.0

import "MethylationBenchmarkPerSampleTasks.wdl" as MethylationBenchmarkPerSampleTasks

workflow MethylationBenchmarkPerSample {
    input {
        File ref
        File fq1gz
        File fq2gz
        File targets
    }

    call MethylationBenchmarkPerSampleTasks.DownsampleReads as DownsampleReads {
        input:
            fq1gz = fq1gz,
            fq2gz = fq2gz
    }

    call MethylationBenchmarkPerSampleTasks.TrimAdapters as TrimAdapters {
        input:
            fq1 = DownsampleReads.fq1Downsampled,
            fq2 = DownsampleReads.fq2Downsampled
    }

    call MethylationBenchmarkPerSampleTasks.FastQC as FastQC {
        input:
            fq1gz = TrimAdapters.fq1Trimmed,
            fq2gz = TrimAdapters.fq2Trimmed
    }

    output {
        File fq1Downsampled = DownsampleReads.fq1Downsampled
        File fq2Downsampled = DownsampleReads.fq2Downsampled
        File fq1Trimmed = TrimAdapters.fq1Trimmed
        File fq2Trimmed = TrimAdapters.fq2Trimmed
        File fastqcReportFq1 = FastQC.htmlReportFq1
        File fastqcReportFq2 = FastQC.htmlReportFq2
    }
}