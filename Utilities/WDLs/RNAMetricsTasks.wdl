version 1.0

task RNASeQC2 {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceAnnotation
    }

    String docker = "gcr.io/broad-cga-aarong-gtex/rnaseqc:latest"
    Int cpu = 1
    Int memory = 64
    Int diskSizeGB = ceil(size(inputBAM, 'GiB') + size(referenceAnnotation, 'GiB')) + 100
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String sampleId = basename(inputBAM, ".bam")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        rnaseqc ~{referenceAnnotation} ~{inputBAM} . --sample ~{sampleId} --mapping-quality 40 --verbose
    >>>

    output {
        File exonCV = "~{sampleId}.exon_cv.tsv"
        File exonReads = "~{sampleId}.exon_reads.gct"
        File geneFragments = "~{sampleId}.gene_fragments.gct"
        File geneReads = "~{sampleId}.gene_reads.gct"
        File geneTPM = "~{sampleId}.gene_tpm.gct"
        File metrics = "~{sampleId}.metrics.tsv"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}
task CollectRNASeqMetrics {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File refFlat
        File ribosomalIntervals
    }

    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    Int cpu = 1
    Int memoryMB = 8000
    Int diskSizeGB = ceil(size(inputBAM, "GiB") + size(referenceGenome, "GiB") + size(referenceGenomeIndex, "GiB")) + 100
    Int javaMemorySize = memoryMB - 1000
    Int maxHeap = memoryMB - 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String outputPrefix = basename(inputBAM, ".bam")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        java -Xms~{javaMemorySize}m -Xmx~{maxHeap}m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
            REF_FLAT=~{refFlat} \
            RIBOSOMAL_INTERVALS=~{ribosomalIntervals} \
            STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
            INPUT=~{inputBAM} \
            VALIDATION_STRINGENCY=SILENT \
            OUTPUT=~{outputPrefix}.rna_metrics
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryMB} MiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

    output {
        File rnaMetrics = outputPrefix + ".rna_metrics"
        File monitoringLog = "monitoring.log"
    }
}