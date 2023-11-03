version 1.0

workflow CrosscheckFingerprints {
    input {
        Array[File] bams
        Array[File] bais
        File haplotypeMap
        String outputPrefix
    }

    scatter(bam in bams) {
        String bamBasenames = basename(bam)
    }

    scatter(bai in bais) {
        String baiBasenames = basename(bai)
    }

    call CrosscheckFingerprintsTask {
        input:
            bams = bams,
            bais = bais,
            haplotypeMap = haplotypeMap,
            outputPrefix = outputPrefix,
            bamBasenames = bamBasenames,
            baiBasenames = baiBasenames
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
        Array[String] bamBasenames
        Array[String] baiBasenames
    }

    Int diskSize = 512
    Int cpu = 4
    Int memory = 64

    command <<<
        for bam in ~{sep=' ' bams}
        do
            mv "${bam}" .
        done;

        for bai in ~{sep=' ' bais}
        do
            mv "${bai}" .
        done;

        java -Xmx60g -jar /usr/picard/picard.jar CrosscheckFingerprints \
            INPUT=~{sep=" INPUT=" bamBasenames} \
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