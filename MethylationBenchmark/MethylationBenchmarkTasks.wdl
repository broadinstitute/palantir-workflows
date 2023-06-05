version 1.0

task DownsampleReads {
    input {
        File fqgz
        Int finalTotalReads = 153930409
        Int rngSeed = 42
        Int cpu = 8
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-dsde-methods/kockan/seqtk@sha256:eb2e9af13f0836fe7652725db4fc82a0e5708778c706bca3fd1f556ecbaba69b"
    }

    String fqBasename = select_first([basename(fqgz, ".fastq.gz"), basename(fqgz, ".fq.gz")])

    command <<<
        gunzip -c ~{fqgz} > ~{fqBasename}.fastq
        seqtk sample -2 -s ~{rngSeed} ~{fqBasename}.fastq ~{finalTotalReads} > ~{fqBasename}.downsampled.fastq
    >>>

    output {
        File fqDownsampled = "~{fqBasename}.downsampled.fastq"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task TrimAdapters {
    input {
        File fq1
        File fq2
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 32
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/trim_galore@sha256:3860476810a6c68c24c504fcaacf0addeca15db3ab207eddf560b86645ae35c5"
    }

    String fq1Basename = basename(fq1, ".fastq.gz")
    String fq2Basename = basename(fq2, ".fastq.gz")

    command <<<
        trim_galore --gzip --cores ~{numThreads} --output_dir . --2colour 20 --paired ~{fq1} ~{fq2}
    >>>

    output {
        File fq1Trimmed = "~{fq1Basename}_val_1.fq.gz"
        File fq2Trimmed = "~{fq2Basename}_val_2.fq.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task FastQC {
    input {
        File fq1
        File fq2
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 32
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/fastqc@sha256:2a3b6bb1df757557cc6daa4a072931b3cd824c1cad4a43c779e70b448a5b1504"
    }

    String fq1Basename = basename(fq1, ".fq.gz")
    String fq2Basename = basename(fq2, ".fq.gz")

    command <<<
        /usr/local/src/FastQC/fastqc --noextract --threads ~{numThreads} --outdir . ~{fq1} ~{fq2}
    >>>

    output {
        File qcFq1 = "~{fq1Basename}_fastqc.html"
        File qcFq2 = "~{fq2Basename}_fastqc.html"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task BWAMethIndex {
    input {
        File ref
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwameth@sha256:20cb5fdf1c1aea1e1209fc0a739d0eec9eef5cb179f5e15887fee83fd7897cc7"
    }

    command <<<
        bwameth.py index ~{ref}

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

task BWAMethAlign {
    input {
        String sampleId
        File ref
        File fq1
        File fq2
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwameth@sha256:20cb5fdf1c1aea1e1209fc0a739d0eec9eef5cb179f5e15887fee83fd7897cc7"
    }

    command <<<
        bwameth.py index ~{ref}

        bwameth.py \
        --reference ~{ref} \
        --threads ~{numThreads} \
        --read-group '@RG\tID:~{sampleId}\tPL:illumina\tSM:~{sampleId}\tLB:~{sampleId}' \
        ~{fq1} ~{fq2} > ~{sampleId}.bam
    >>>

    output {
        File bam = "~{sampleId}.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SAMBamba {
    input {
        File ref
        File bam
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-dsde-methods/kockan/sambamba@sha256:a27ab0121ffb3b5a5346ddb0d531a90966923015e8a945de26d2465f3103da73"
    }

    String bamBasename = basename(bam, ".bam")
    #String sortedBAM = bamBasename + ".sorted.bam"

    command <<<
        sambamba view -h \
            -t ~{numThreads} \
            -T ~{ref} \
            -F 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -o ~{bamBasename}.filtered.bam \
            ~{bam}
    >>>

    output {
        #File sortedBam = "~{sortedBAM}"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

#        sambamba sort \
#            -t ~{numThreads} \
#            -m ~{memoryGB}GiB \
#            --tmpdir . \
#            -o /dev/stdout \
#            temp.bam | sambamba view -h -t ~{numThreads} -o ~{sortedBAM} -T ~{ref} -f bam /dev/stdin

task IndexBAM {
    input {
        File bam
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 64
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-dsde-methods/kockan/samtools@sha256:b0f4520282c18967e279071615dcc7685ee9457649928664d68728add6f01156"
    }

    command <<<
        samtools index -@ ~{numThreads} ~{bam} > ~{bam}.bai
    >>>

    output {
        File bai = "~{bam}.bai"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task MarkDuplicates {
    input {
        String sampleId
        File ref
        File bam
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    String bamBasename = basename(bam, ".bam")

    command <<<
        java -Xmx4g -Xms4g -jar /usr/picard/picard.jar MarkDuplicates \
            -I ~{bam} \
            -O ~{bamBasename}.markdup.bam \
            -R ~{ref} \
            -M ~{sampleId}.sorted.picard_markdup_raw_metrics \
            --CREATE_INDEX false \
            --MAX_RECORDS_IN_RAM 1000 \
            --SORTING_COLLECTION_SIZE_RATIO 0.15 \
            --ASSUME_SORT_ORDER coordinate \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
    >>>

    output {
        File outputBam = "~{bamBasename}.markdup.bam"
        File outputMetrics = "~{sampleId}.sorted.picard_markdup_raw_metrics"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task ConvertBedToIntervalList {
    input {
        File bed
        File dict
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    command <<<
        java -Xmx4g -Xms4g -jar /usr/picard/picard.jar BedToIntervalList \
            -I ~{bed} \
            -O ~{bed}.intervals \
            -SD ~{dict}
    >>>

    output {
        File intervalList = "~{bed}.intervals"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task CollectHsMetrics {
    input {
        File ref
        File bam
        File intervals
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    command <<<
        java -Xmx4g -Xms4g -jar /usr/picard/picard.jar CollectHsMetrics \
            -I ~{bam} \
            -O metrics.txt \
            -R ~{ref} \
            --BAIT_INTERVALS ~{intervals} \
            --TARGET_INTERVALS ~{intervals} \
            --MINIMUM_MAPPING_QUALITY 20 \
            --COVERAGE_CAP 1000 \
            --PER_TARGET_COVERAGE asd.text \
            --NEAR_DISTANCE 500
    >>>

    output {
        File metrics = "metrics.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task CollectMultipleMetrics {
    input {
        File ref
        File bam
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:61880f4b6190a955d30ef61b3e3d48b1889974765dea64ee793450bf97144389"
    }

    command <<<
        java -Xmx4g -Xms4g -jar /usr/picard/picard.jar CollectMultipleMetrics \
            -I ~{bam} \
            -O rawmetrics.txt \
            -R ~{ref} \
            --PROGRAM null \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics
    >>>

    output {
        File rawMetrics = "rawmetrics.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task CallMethylation {
    input {
        String sampleId
        File bam
        File ref
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/methyldackel@sha256:36a7349df6bad066df5af7db33b913fe902d6ae3a67bb1544552b7c8a2da90a5"
    }

    command <<<
        MethylDackel mbias ~{ref} ~{bam} ~{sampleId}

        MethylDackel extract \
        --minDepth 10 \
        --maxVariantFrac 0.25 \
        --OT X,X,X,X \
        --OB X,X,X,X \
        --mergeContext ~{ref} ~{bam} \
        -o ~{sampleId}

        MethylDackel extract \
        --minDepth 10 \
        --maxVariantFrac 0.25 \
        --OT 0,0,0,98 \
        --OB 0,0,3,0 \
        --cytosine_report \
        --CHH --CHG ~{ref} ~{bam} \
        -o ~{sampleId}_report
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

#task CollectMethylationStatistics {
#    input {
#        File bam
#        Int cpu = 4
#        Int numThreads = 8
#        Int memoryGB = 32
#        Int diskSizeGB = 256
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/samtools@sha256:b0f4520282c18967e279071615dcc7685ee9457649928664d68728add6f01156"
#    }
#
#    command <<<
#        samtools stats ~{bam} | grep ^SN | cut -f 2-
#    >>>
#
#    output {
#    }
#
#    runtime {
#        cpu: cpu
#        memory: "~{memoryGB} GiB"
#        disks: "local-disk ~{diskSizeGB} HDD"
#        docker: docker
#    }
#}