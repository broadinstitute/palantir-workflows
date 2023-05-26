version 1.0

task FlamesTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/flames@sha256:e9b5d5152179e1a820afde3b147586a8ce7440738bf456af74b22ca4cfa0e8cb"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir fq

        samtools fastq ~{inputBAM} > ./fq/temp.fastq

        python3 /usr/local/src/FLAMES/python/bulk_long_pipeline.py \
        --gff3 ~{referenceAnnotation} \
        --genomefa ~{referenceGenome} \
        --fq_dir ./fq \
        --inbam ~{inputBAM} \
        --outdir .
    >>>

    output {
        File flamesGFF = "isoform_annotated.gff3"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow Flames {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
    }

    call FlamesTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    output {
        File flamesGFF = FlamesTask.flamesGFF
        File monitoringLog = FlamesTask.monitoringLog
    }
}