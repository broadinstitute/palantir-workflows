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

    #call MethylationBenchmarkPerSampleTasks.FastQC as FastQC {
    #    input:
    #        fq1 = fq1,
    #        fq2 = fq2
    #}
    
    output {
        File fq1Downsampled = DownsampleReads.fq1Downsampled
        File fq2Downsampled = DownsampleReads.fq2Downsampled
    }
}