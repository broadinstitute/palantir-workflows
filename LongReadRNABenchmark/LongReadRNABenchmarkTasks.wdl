version 1.0

task IsoQuant {
    input {
        File inputBAM
        File referenceGenome
        File reducedAnnotation
        String datasetName
        String dataType
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 250

    command {
        /usr/local/src/IsoQuant-3.1.1/isoquant.py \
        --reference ~{referenceGenome} \
        --complete_genedb \
        --genedb ~{reducedAnnotation} \
        --bam ~{inputBAM} \
        -d ~{dataType} \
        -t ~{numThreads} \
        -l ~{datasetName} \
        -o "IsoQuant_out_~{datasetName}"
    }

    output {
        File isoQuantGTF = "IsoQuant_out_~{datasetName}.transcript_models.gtf"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}