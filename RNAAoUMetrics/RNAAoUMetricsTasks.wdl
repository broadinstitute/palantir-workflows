version 1.0

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

    String outputPrefix = basename(alignment, ".bam")

    command <<<
        java -Xmx32g -Xms4g -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
        --INPUT ~{alignment} \
        --OUTPUT ~{outputPrefix}.rna_metrics
        --REF_FLAT ~{refFlat} \
        --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
        --RIBOSOMAL_INTERVALS ~{ribosomalIntervals} \
        --VALIDATION_STRINGENCY SILENT \
    >>>

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

    output {
        File rnaSeqMetrics = outputPrefix + ".rna_metrics"
    }
}