version 1.0

task CollectInsertSizeMetrics {
    input {
        File alignment
        File alignmentIndex
        File reference
        File referenceIndex
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    }

    String outputPrefix = basename(alignment, ".cram")

    command <<<
        java -Xmx32g -Xms4g -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
        --INPUT ~{alignment} \
        --REFERENCE_SEQUENCE ~{reference} \
        --OUTPUT ~{outputPrefix}.insert_size_metrics.txt \
        --Histogram_FILE ~{outputPrefix}.insert_size_histogram.pdf \
        --VALIDATION_STRINGENCY SILENT
    >>>

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

    output {
        File insertSizeMetrics = outputPrefix + ".insert_size_metrics.txt"
        File insertSizeHistogram = outputPrefix + ".insert_size_histogram.pdf"
    }
}

task CollectAlignmentSummaryMetrics {
    input {
        File alignment
        File alignmentIndex
        File reference
        File referenceIndex
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    }

    String outputPrefix = basename(alignment, ".cram")

    command <<<
        java -Xmx32g -Xms4g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
        --INPUT ~{alignment} \
        --REFERENCE_SEQUENCE ~{reference} \
        --OUTPUT ~{outputPrefix}.alignment_summary_metrics.txt \
        --VALIDATION_STRINGENCY SILENT
    >>>

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

    output {
        File alignmentSummaryMetrics = outputPrefix + ".alignment_summary_metrics.txt"
    }
}

task CollectRNASeqMetrics {
    input {
        File alignment
        File alignmentIndex
        File reference
        File referenceIndex
        File refFlat
        File ribosomalIntervals
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    }

    String outputPrefix = basename(alignment, ".cram")

    command <<<
        java -Xmx32g -Xms4g -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
        --INPUT ~{alignment} \
        --REFERENCE_SEQUENCE ~{reference} \
        --OUTPUT ~{outputPrefix}.rnaseq_metrics.txt \
        --REF_FLAT ~{refFlat} \
        --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
        --RIBOSOMAL_INTERVALS ~{ribosomalIntervals} \
        --VALIDATION_STRINGENCY SILENT
    >>>

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

    output {
        File rnaSeqMetrics = outputPrefix + ".rnaseq_metrics.txt"
    }
}

task MarkDuplicates {
    input {
        File alignment
        File alignmentIndex
        File reference
        File referenceIndex
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 256
        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    }

    String outputPrefix = basename(alignment, ".cram")

    command <<<
        java -Xmx32g -Xms4g -jar /usr/picard/picard.jar MarkDuplicates \
        --INPUT ~{alignment} \
        --REFERENCE_SEQUENCE ~{reference} \
        --METRICS_FILE ~{outputPrefix}.markdup_metrics.txt \
        --OUTPUT /dev/null \
        --VALIDATION_STRINGENCY SILENT
    >>>

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

    output {
        File duplicationMetrics = outputPrefix + ".markdup_metrics.txt"
    }
}