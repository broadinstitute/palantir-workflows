version 1.0

task DownsampleReads {
    input {
        File fq1gz
        File fq2gz
        Int finalTotalReads = 153930409
        Int rngSeed = 42
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/seqtk@sha256:eb2e9af13f0836fe7652725db4fc82a0e5708778c706bca3fd1f556ecbaba69b"
    }

    String fq1Basename = select_first([basename(fq1gz, ".fastq.gz"), basename(fq1gz, ".fq.gz")])
    String fq2Basename = select_first([basename(fq2gz, ".fastq.gz"), basename(fq2gz, ".fq.gz")])

    command <<<
        gunzip -c ~{fq1gz} > ~{fq1Basename}.fastq
        gunzip -c ~{fq2gz} > ~{fq2Basename}.fastq

        seqtk sample -2 -s ~{rngSeed} ~{fq1Basename}.fastq ~{finalTotalReads} > ~{fq1Basename}.downsampled.fastq
        seqtk sample -2 -s ~{rngSeed} ~{fq2Basename}.fastq ~{finalTotalReads} > ~{fq2Basename}.downsampled.fastq
    >>>

    output {
        File fq1Downsampled = "~{fq1Basename}.downsampled.fastq"
        File fq2Downsampled = "~{fq2Basename}.downsampled.fastq"
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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/trim_galore@sha256:3860476810a6c68c24c504fcaacf0addeca15db3ab207eddf560b86645ae35c5"
    }

    String fq1Basename = basename(fq1, ".downsampled.fastq")
    String fq2Basename = basename(fq2, ".downsampled.fastq")

    command <<<
        trim_galore --gzip --cores ~{cpu} --output_dir . --2colour 20 --paired ~{fq1} ~{fq2}
    >>>

    output {
        File fq1Trimmed = "~{fq1Basename}.trimmed.fq.gz"
        File fq2Trimmed = "~{fq2Basename}.trimmed.fq.gz"
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
        File fq1gz
        File fq2gz
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/fastqc@sha256:ec113d537e232de7b8030a9944f16274ca855f487ecdd8e2f8db1aebecfeeeb5"
    }

    String fq1Basename = select_first([basename(fq1gz, ".fastq.gz"), basename(fq1gz, ".fq.gz")])
    String fq2Basename = select_first([basename(fq2gz, ".fastq.gz"), basename(fq2gz, ".fq.gz")])

    command <<<
        fastqc --noextract --threads ~{numThreads} ~{fq1gz} ~{fq2gz}
    >>>

    output {
        File htmlReportFq1 = "~{fq1Basename}_fastqc.html"
        File htmlReportFq2 = "~{fq2Basename}_fastqc.html"
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
        File amb
        File ann
        File bwt
        File pac
        File sa
        File fq1
        File fq2
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bwameth@sha256:c2415b900b8aa96d39c0c9e615e6f6eece37bf142bdbfee1ec783434390c34d9"
    }

    command <<<
        cp ~{amb} ~{ann} ~{bwt} ~{pac} ~{sa} .

        bwameth.py \
            --reference ~{ref} \
            --threads ~{numThreads} \
            --read-group '@RG\tSAMPLE_ID:1\tPL:illumina\tLB:SAMPLE\tSM:SAMPLE' \
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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/sambamba@sha256:a27ab0121ffb3b5a5346ddb0d531a90966923015e8a945de26d2465f3103da73"
    }

    String bamBasename = basename(bam, ".bam")
    String sortedBAM = bamBasename + ".sorted.bam"

    command <<<
        sambamba view -h \
            -t 16 \
            -T ~{ref} \
            --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam \
            -l 0 ~{bam} \
            -o temp.bam

        sambamba sort \
            -t 16 \
            -m 30Gib \
            --tmpdir /tmp/ \
            -o /dev/stdout \
            -l temp.bam | sambamba view -h -t 16 -o ~{sortedBAM} -T ~{ref} -f bam /dev/stdin
    >>>

    output {
        File sortedBam = "~{sortedBAM}"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

#task IndexBAM {
#    input {
#        File bam
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        samtools index -@ 16 ~{bam} > ~{bam}.bai
#    >>>
#
#    output {
#        File bai = "~{bam}.bai"
#    }
#
#    runtime {
#        cpu: cpu
#        memory: "~{memoryGB} GiB"
#        disks: "local-disk ~{diskSizeGB} HDD"
#        docker: docker
#    }
#}
#
#task MarkDuplicates {
#    input {
#        File reference
#        File bam
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    String bamBasename = basename(bam, ".bam")
#
#    command <<<
#        java -Xmx4g -Xms4g -jar picard.jar MarkDuplicates \
#        -I ~{bam} \
#        -O ~{bamBasename}.markdup.bam \
#        -R ~{reference} \
#        -M S1_methylome_sorted.picard_markdup_raw_metrics \ #TODO: properly construct this name for the sample!
#        --CREATE_INDEX false \
#        --MAX_RECORDS_IN_RAM 1000 \
#        --SORTING_COLLECTION_SIZE_RATIO 0.15 \
#        --ASSUME_SORT_ORDER coordinate \
#        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
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
#
#task ConvertBedToIntervalList {
#    input {
#        File bed
#        File dict
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        java -Xmx4g -Xms4g -jar picard.jar BedToIntervalList \
#        -I ~{bed} \
#        -O ~{bed}.intervals \
#        SD ~{dict}
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
#
#task CollectHsMetrics {
#    input {
#        File bam
#        File intervals
#        File reference
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        java -Xmx4g -Xms4g -jar picard.jar CollectHsMetrics \
#        -I ~{bam} \
#        -O metrics.txt \
#        -R ~{reference} \
#        --BAIT_INTERVALS ~{intervals} \
#        --TARGET_INTERVALS ~{intervals} \
#        --MINIMUM_MAPPING_QUALITY 20 \
#        --COVERAGE_CAP 1000 \
#        --PER_TARGET_COVERAGE asd.text \
#        --NEAR_DISTANCE 500
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
#
#task CollectMultipleMetrics {
#    input {
#        File bam
#        File reference
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        java -Xmx4g -Xms4g -jar picard.jar CollectMultipleMetrics \
#        -I ~{bam} \
#        -O rawmetrics.txt \
#        -R ~{reference} \
#        --PROGRAM null \
#        --PROGRAM CollectGcBiasMetrics \
#        --PROGRAM CollectInsertSizeMetrics \
#        --PROGRAM CollectAlignmentSummaryMetrics
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
#
#task CallMethylation {
#    input {
#        String sample
#        File bam
#        File reference
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        MethylDackel mbias ~{reference} ~{bam} ~{sample}
#
#        MethylDackel extract \
#        --minDepth 10 \
#        --maxVariantFrac 0.25 \
#        --OT X,X,X,X \
#        --OB X,X,X,X \
#        --mergeContext ~{reference} ~{bam} \
#        -o ~{sample}
#
#        MethylDackel extract \
#        --minDepth 10 \
#        --maxVariantFrac 0.25 \
#        --OT 0,0,0,98 \
#        --OB 0,0,3,0 \
#        --cytosine_report \
#        --CHH --CHG ~{reference} ~{bam} \
#        -o ~{sample}_report
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
#
#task CollectMethylationStatistics {
#    input {
#        File bam
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
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