version 1.0

task IsoQuantTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation
        String datasetName
        String dataType
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant@sha256:9bd8cd8c3a04e02599e10e0b484127fb763a39499302d4c859d230942f9a2d15"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = if defined(referenceAnnotation) then "IsoQuant_out_~{datasetName}" else "IsoQuant_denovo_out_~{datasetName}"
    String completeGeneDBOption = if defined(referenceAnnotation) then "--complete_genedb" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.2.0/isoquant.py \
        --reference ~{referenceGenome} \
        ~{completeGeneDBOption} \
        ~{"--genedb " + referenceAnnotation} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output ~{outputPrefix}
    >>>

    output {
        File isoQuantGTF = "~{outputPrefix}/~{datasetName}/~{datasetName}.transcript_models.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow IsoQuant {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation
        String datasetName
        String dataType
    }

    call IsoQuantTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
            dataType = dataType
    }

    output {
        File isoQuantGTF = IsoQuantTask.isoQuantGTF
        File monitoringLog = IsoQuantTask.monitoringLog
    }
}