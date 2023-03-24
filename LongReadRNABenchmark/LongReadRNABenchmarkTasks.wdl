version 1.0

task IsoQuant {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        String dataType
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String referenceAnnotationBasename = basename(referenceAnnotation, ".reduced.gtf")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.1.1/isoquant.py \
        --reference ~{referenceGenome} \
        --complete_genedb \
        --genedb ~{referenceAnnotation} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output "IsoQuant_out_~{datasetName}"
    >>>

    output {
        File isoQuantGTF = "IsoQuant_out_~{datasetName}/~{datasetName}/~{datasetName}.transcript_models.gtf"
        File isoQuantDB = "IsoQuant_out_~{datasetName}/~{referenceAnnotationBasename}.reduced.db"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task StringTie {
    input {
        File inputBAM
        File referenceAnnotation
        String datasetName
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/stringtie \
        -o "StringTie_out_~{datasetName}.gtf" \
        -G ~{referenceAnnotation} \
        -p ~{numThreads} \
        -L \
        ~{inputBAM} \
    >>>

    output {
        File stringTieGTF = "StringTie_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReducedAnnotationGFFCompare {
    input {
        File reducedAnnotationDB
        File isoQuantGTF
        File stringTieGTF
        String datasetName
    }

    # TODO: Need to add gffcompare to the docker images
    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 300
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String reducedAnnotationPrefix = basename(reducedAnnotationDB, ".reduced.db")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{isoQuantGTF} \
        --tool "isoquant" \
        --output "~{datasetName}_isoquant_reduced_db"

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{stringTieGTF} \
        --tool "stringtie" \
        --output "~{datasetName}_stringtie_reduced_db"
    >>>

    output {
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}