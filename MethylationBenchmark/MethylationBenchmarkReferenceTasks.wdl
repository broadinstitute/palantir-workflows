version 1.0

task GenerateBWAMethIndex {
    input {
        File ref
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwameth@sha256:20cb5fdf1c1aea1e1209fc0a739d0eec9eef5cb179f5e15887fee83fd7897cc7"
    }

    String refBasename = basename(ref)
    String c2t = "~{refBasename}.bwameth.c2t"
    String amb = "~{refBasename}.bwameth.c2t.amb"
    String ann = "~{refBasename}.bwameth.c2t.ann"
    String bwt = "~{refBasename}.bwameth.c2t.bwt"
    String pac = "~{refBasename}.bwameth.c2t.pac"
    String sa = "~{refBasename}.bwameth.c2t.sa"

    command <<<
        mv ~{ref} .
        bwameth.py index ~{refBasename}
        mkdir ~{refBasename}.bwameth_index
        mv ~{c2t} ~{amb} ~{ann} ~{bwt} ~{pac} ~{sa} ~{refBasename}.bwameth_index
        tar -zcvf ~{refBasename}.bwameth_index.tar.gz ~{refBasename}.bwameth_index
    >>>

    output {
        File bwamethIndex = "~{refBasename}.bwameth_index.tar.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GenerateFASTAIndex {
    input {
        File ref
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us.gcr.io/broad-dsde-methods/kockan/samtools@sha256:b0f4520282c18967e279071615dcc7685ee9457649928664d68728add6f01156"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} .
        samtools faidx ~{refBasename} > ~{refBasename}.fai
    >>>

    output {
        File fai = "~{refBasename}.fai"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task CreateSequenceDictionary {
    input {
        File ref
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    String refBasename = basename(ref)

    command <<<
        java -Xmx32g -Xms8g -jar /usr/picard/picard.jar CreateSequenceDictionary \
        REFERENCE=~{ref} \
        OUTPUT=~{refBasename}.dict
    >>>

    output {
        File dict = "~{refBasename}.dict"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}