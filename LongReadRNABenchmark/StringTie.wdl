version 1.0

task StringTieTask {
    input {
        File inputBAM
        File? referenceAnnotation
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie@sha256:ca2a163c7acdcacba741ea98d573080c15f153de18bd1566d24e8d2f1729ce89"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = if defined(referenceAnnotation) then "StringTie_out_~{datasetName}" else "StringTie_denovo_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        stringtie \
        -o "~{outputPrefix}.gtf" \
        ~{"-G " + referenceAnnotation} \
        -p ~{numThreads} \
        -L ~{inputBAM}
    >>>

    output {
        File stringTieGTF = "~{outputPrefix}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow StringTie {
    input {
        File inputBAM
        File? referenceAnnotation
        String datasetName
    }

    call StringTieTask {
        input:
            inputBAM = inputBAM,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    output {
        File stringTieGTF = StringTieTask.stringTieGTF
        File monitoringLog = StringTieTask.monitoringLog
    }
}