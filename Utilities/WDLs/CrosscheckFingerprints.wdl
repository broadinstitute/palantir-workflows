version 1.0

workflow CrosscheckFingerprints {
    input {
        Array[File] bams
        Array[File] bais
        File haplotypeMap
        String outputPrefix
    }

    call CrosscheckFingerprintsTask {
        input:
            bams = bams,
            bais = bais,
            haplotypeMap = haplotypeMap,
            outputPrefix = outputPrefix
    }

    output {
        File outputMetrics = CrosscheckFingerprintsTask.outputMetrics
    }
}

task CrosscheckFingerprintsTask {
    input {
        Array[File] bams
        Array[File] bais
        File haplotypeMap
        String outputPrefix
    }

    Int diskSize = 512
    Int cpu = 4
    Int memory = 64

    command <<<
        java -Xmx63g -jar /usr/picard/picard.jar CrosscheckFingerprints \
            INPUT=~{sep=" INPUT=" bams} \
            HAPLOTYPE_MAP=~{haplotypeMap} \
            CROSSCHECK_BY=SAMPLE \
            OUTPUT=~{outputPrefix}.crosscheck_metrics
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:3.0.0"
        disks: "local-disk " + diskSize + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File outputMetrics = "~{outputPrefix}.crosscheck_metrics"
    }
}