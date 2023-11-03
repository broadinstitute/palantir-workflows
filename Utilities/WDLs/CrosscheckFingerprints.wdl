version 1.0

workflow CrosscheckFingerprints {
    input {
        Array[String] bams
        Array[String] bais
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
        Array[String] bams
        Array[String] bais
        File haplotypeMap
        String outputPrefix
    }

    Int diskSize = 512
    Int cpu = 4
    Int memory = 64

    command <<<
        java -jar /usr/picard/picard.jar CrosscheckFingerprints \
            INPUT=~{write_lines(bams)} \
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