version 1.0

struct MetricsFiles {
    File error_by_gc
    File error_by_read_ordinality
    File error_by_ref_base
    File error_by_pre_dinuc
    File error_by_post_dinuc
    File error_by_cycle
    File error_by_insert_length
    File error_by_base_quality
    File error_by_mapping_quality
    File error_by_one_base_padded_context
    File error_by_indel_length
}

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
        MetricsFiles error_metrics = CollectErrorMetrics.error_metrics
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
        gatk CollectSamErrorMetrics -I ~{bam} -V ~{vcf} -R ~{reference_fasta} -O ~{output_basename} \
            --ERROR_METRICS ERROR:REFERENCE_BASE \
            --ERROR_METRICS ERROR:PRE_DINUC \
            --ERROR_METRICS ERROR:POST_DINUC
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        disks: "local-disk " + disk_size + " HDD"
        cpu: 4
        memory: mem_gb + " GB"
    }

    output {
        MetricsFiles error_metrics = {
                "error_by_gc": output_basename + ".error_by_gc",
                "error_by_read_ordinality": output_basename + ".error_by_read_ordinality",
                "error_by_ref_base": output_basename + ".error_by_ref_base",
                "error_by_pre_dinuc": output_basename + ".error_by_pre_dinuc",
                "error_by_post_dinuc": output_basename + ".error_by_post_dinuc",
                "error_by_cycle": output_basename + ".error_by_cycle",
                "error_by_insert_length": output_basename + ".error_by_insert_length",
                "error_by_base_quality": output_basename + ".error_by_base_quality",
                "error_by_mapping_quality": output_basename + ".error_by_mapping_quality",
                "error_by_one_base_padded_context": output_basename + ".error_by_one_base_padded_context",
                "error_by_indel_length": output_basename + ".error_by_indel_length",
            }
    }
}