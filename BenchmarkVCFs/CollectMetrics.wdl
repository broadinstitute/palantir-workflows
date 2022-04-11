version 1.0

workflow CollectMetrics {
    input {
        File vcf
        File vcf_index
        File bam
        File bam_index

        File reference_fasta
        File reference_index
    }

    call CollectErrorMetrics {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            bam = bam,
            bam_index = bam_index,
            reference_fasta = reference_fasta,
            reference_index = reference_index
    }

    output {
        Array[File] error_metrics = CollectErrorMetrics.error_metrics
    }
}

task CollectErrorMetrics {
    input {
        File vcf
        File vcf_index
        File bam
        File bam_index

        File reference_fasta
        File reference_index

        Int preemptible = 1
        String docker = "us.gcr.io/broad-gatk/gatk:4.2.6.0"
        Int mem_gb = 4
    }

    String output_basename = sub(basename(bam), "\.(bam|cram)$", "")

    Int disk_size = ceil(1.5 * (size(bam,"GB") + size(bam_index,"GB") + size(vcf_index,"GB") + size(vcf_index,"GB")) + size(reference_fasta, "GB") + size(reference_index, "GB"))

    command <<<
        gatk CollectSamErrorMetrics -I ~{bam} -V ~{vcf} -R ~{reference_fasta} -O output_basename \
            --ERROR_METRICS \
                ERROR:GC_CONTENT \
                ERROR:READ_ORDINALITY \
                ERROR:REFERENCE_BASE \
                ERROR:PRE_DINUC \
                ERROR:POST_DINUC \
                ERROR:CYCLE \
                ERROR:INSERT_LENGTH \
                ERROR:BASE_QUALITY \
                ERROR:MAPPING_QUALITY \
                ERROR:ONE_BASE_PADDED_CONTEXT \
                ERROR:INDEL_LENGTH
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        disks: "local-disk " + disk_size + " HDD"
        cpu: 4
        memory: mem_gb + " GB"
    }

    output {
        Array[File] error_metrics = glob(output_basename + ".*")
    }
}