version 1.0

task CreateBigWigFromBedGraph {
    input {
        String sampleId
        File bedGraph
        File chromSizes
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bedgraphtobigwig@sha256:4d57ab760cf6f79df1f07b7cfb12c93afd4dcff9ae20dfe767d16dd958a157e9"
    }

    command <<<
        awk '{print $1,$2,$3,$4}' ~{bedGraph} | sort -k1,1 -k2,2n | grep -v "track" > temp.bedGraph
        bedGraphToBigWig temp.bedGraph ~{chromSizes} ~{sampleId}.bw
    >>>

    output {
        File bigWig = "~{sampleId}.bw"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task MultiBigWigSummary {
    input {
        Array[String] sampleIds
        Array[File] bigWigs
        String outputPrefix
        Int binSize = 1000
        Int cpu = 1
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