version 1.0

task MultiBigWigSummary {
    input {
        Array[String] sampleIds
        Array[File] bigWigs
        String outputPrefix
        Int binSize = 1000
        Int cpu = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us.gcr.io/broad-dsde-methods/kockan/deeptools@sha256:79a66b54fb00828c24208e8f365214198d4f7b4696a69d8baad0e3b0daf8c412"
    }

    command <<<
        multiBigwigSummary bins \
        --binSize 1000 \
        -b ~{sep=" " bigWigs} \
        --labels ~{sep=" " sampleIds} \
        -o ~{outputPrefix}.npz
    >>>

    output {
        File multiBigWigSummary = "~{outputPrefix}.npz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}
    
workflow MethylationSupplementarySet {
    input {
        Array[String] sampleIds
        Array[File] bigWigs
        String outputPrefix
        Int binSize
    }

    call MultiBigWigSummary {
        input:
            sampleIds = sampleIds,
            bigWigs = bigWigs,
            outputPrefix = outputPrefix,
            binSize = binSize
    }

    output {
        File multiBigWigSummary = MultiBigWigSummary.multiBigWigSummary
    }
}

