version 1.0

import "MethylationBenchmarkPerSampleTasks.wdl" as MethylationBenchmarkPerSampleTasks

workflow MethylationBenchmark {
    input {
        File ref
        File fq1gz
        File fq2gz
        File targets
    }

    #call MethylationBenchmarkTasks.FastQC as FastQC {
    #    input:
    #        fq1 = fq1,
    #        fq2 = fq2
    #}

    output {
    }
}