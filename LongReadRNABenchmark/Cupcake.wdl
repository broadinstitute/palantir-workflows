version 1.0

task CupcakeTask {
    input {
        File inputBAM
        File inputBAMIndex
        String datasetName
        Int cpu = 16
        Int memoryGB = 512
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/cdna-cupcake@sha256:fca085dde170c995b5691d07aae2d56ab4426b7651a913957f029e628a0167c2"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = "Cupcake_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        python3 /usr/local/src/remove_fastq_duplicates.py temp.fastq

        python3 /usr/local/src/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py \
        --input out.fastq --fq \
        --bam ~{inputBAM} \
        --prefix ~{outputPrefix} \
        --cpus 1
    >>>

    output {
        File cupcakeGFF = "~{outputPrefix}.collapsed.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow Cupcake {
    input {
        File inputBAM
        File inputBAMIndex
        String datasetName

    }

    call CupcakeTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            datasetName = datasetName
    }

    output {
        File cupcakeGFF = CupcakeTask.cupcakeGFF
        File monitoringLog = CupcakeTask.monitoringLog
    }
}