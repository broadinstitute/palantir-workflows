version 1.0

task IsoSeqTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoseq3@sha256:e715dda61f295d6825c0f4bea5133d583158db5d63c550ff186ee59f1ff10385"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        pbmm2 align --num-threads ~{numThreads} --preset ISOSEQ --sort ~{referenceGenome} temp.fastq  pbmm_realigned.bam

        isoseq3 collapse pbmm_realigned.bam "IsoSeq_out_~{datasetName}.gff"
    >>>

    output {
        File isoSeqGFF = "IsoSeq_out_~{datasetName}.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow IsoSeq {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        String datasetName
    }

    call IsoSeqTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            datasetName = datasetName
    }

    output {
        File isoSeqGFF = IsoSeqTask.isoSeqGFF
        File monitoringLog = IsoSeqTask.monitoringLog
    }
}