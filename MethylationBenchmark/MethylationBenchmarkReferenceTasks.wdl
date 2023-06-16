version 1.0

task GenerateBWAMethIndex {
    input {
        File ref
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwameth@sha256:20cb5fdf1c1aea1e1209fc0a739d0eec9eef5cb179f5e15887fee83fd7897cc7"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} .

        ls -lha

        bwameth.py index ~{refBasename}

        ls -lha
    >>>

    output {
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GenerateBWAIndex {
    input {
        File ref
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwa@sha256:2b87ef641b91fbb1c45f852d09216dec4122a1ca5b525ee2dea3d87330b0004f"
    }

    command <<<
        bwa index -a bwtsw ~{ref}
    >>>

    output {
        File amb = "~{ref}.amb"
        File ann = "~{ref}.ann"
        File bwt = "~{ref}.bwt"
        File pac = "~{ref}.pac"
        File sa = "~{ref}.sa"
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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/samtools@sha256:b0f4520282c18967e279071615dcc7685ee9457649928664d68728add6f01156"
    }

    command <<<
        samtools faidx ~{ref} > ~{ref}.fai
    >>>

    output {
        File fai = "~{ref}.fai"
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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    String refBasename = basename(ref, ".fa")

    command <<<
        java -Xmx4g -Xms4g -jar /usr/picard/picard.jar CreateSequenceDictionary \
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