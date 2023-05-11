version 1.0

#task FastQC {
#    input {
#        File fq1
#        File fq2
#        Int cpu = 8
#        Int numThreads = 16
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        fastqc --noextract --threads ~{numThreads} ~{fq1} ~{fq2}
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
#task TrimAdapters {
#    input {
#        File fq1
#        File fq2
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    String fq1Basename = basename(fq1, "_downsampled.fastq")
#    String fq2Basename = basename(fq2, "_downsampled.fastq")
#
#    command <<<
#        trim_galore --gzip --cores ~{cpu} --output_dir . --2colour 20 --paired ~{fq1} ~{fq2}
#    >>>
#
#    output {
#        File fq1Trimmed = "~{fq1Basename}.trimmed.fq.gz"
#        File fq2Trimmed = "~{fq2Basename}.trimmed.fq.gz"
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
#task BWAMeth {
#    input {
#        File reference
#        File fq1
#        File fq2
#        Int cpu = 16
#        Int numThreads = 32
#        Int memoryGB = 64
#        Int diskSizeGB = 100
#        String docker = "us.gcr.io/broad-dsde-methods/kockan/..."
#    }
#
#    command <<<
#        bwameth.py --reference ~{reference} -t 16 --read-group '@RG\tSAMPLE_ID:1\tPL:illumina\tLB:SAMPLE\tSM:SAMPLE' ~{fq1} ~{fq2}
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
#task SAMBamba {
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
#    String sortedBAM = bamBasename + "_sorted.bam"
#
#    command <<<
#        sambamba view -h \
#        -t 16 \
#        -T ~{reference} \
#        --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
#        -f bam \
#        -l 0 ~{bam} \
#        -o /tmp/tmpzs73858y.bam
#
#        sambamba sort -t 16 -m 30Gib --tmpdir /tmp/ -o /dev/stdout -l 0 /tmp/tmpzs73858y.bam | sambamba view -h -t 16 -o ~{sortedBAM} -T ~{reference} -f bam /dev/stdin
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