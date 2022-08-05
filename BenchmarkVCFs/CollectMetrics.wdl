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
}

workflow CollectMetrics {
    input {
        File vcf
        File vcf_index
        File bam
        File bam_index

        File intervals_for_wgs_metrics
        File? intervals_for_error_metrics
        Int read_length

        Boolean use_fast_algorithm = false

        File reference_fasta
        File reference_index
        File reference_dict

        String gatk_tag
        Int preemptible = 1
    }

    call CollectErrorMetrics {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            bam = bam,
            bam_index = bam_index,
            intervals = intervals_for_error_metrics,
            reference_fasta = reference_fasta,
            reference_index = reference_index,
            reference_dict = reference_dict,
            gatk_tag = gatk_tag,
            preemptible = preemptible
    }

    call CollectWgsMetrics {
        input:
            input_bam = bam,
            input_bam_index = bam_index,
            wgs_coverage_interval_list = intervals_for_wgs_metrics,
            read_length = read_length,
            use_fast_algorithm = use_fast_algorithm,
            reference_fasta = reference_fasta,
            reference_index = reference_index,
            reference_dict = reference_dict,
            gatk_tag = gatk_tag,
            preemptible = preemptible
    }

    output {
        MetricsFiles error_metrics = CollectErrorMetrics.error_metrics
        File wgs_metrics = CollectWgsMetrics.metrics
    }
}

task CollectErrorMetrics {
    input {
        File vcf
        File vcf_index
        File bam
        File bam_index

        File? intervals

        File reference_fasta
        File reference_index
        File reference_dict

        Int preemptible
        String gatk_tag
        Int mem_gb = 4
    }

    String output_basename = sub(basename(bam), "\.(bam|cram)$", "")

    Int disk_size = ceil(1.5 * (size(bam,"GB") + size(bam_index,"GB") + size(vcf_index,"GB") + size(vcf_index,"GB")) + size(reference_fasta, "GB") + size(reference_index, "GB"))

    command <<<
        gatk CollectSamErrorMetrics -I ~{bam} -V ~{vcf} -R ~{reference_fasta} -O ~{output_basename} \
            ~{true="-L" false="" defined(intervals)} ~{default="" intervals} \
            --ERROR_METRICS null \
            --ERROR_METRICS ERROR:GC_CONTENT \
            --ERROR_METRICS ERROR:READ_ORDINALITY \
            --ERROR_METRICS ERROR:REFERENCE_BASE \
            --ERROR_METRICS ERROR:PRE_DINUC \
            --ERROR_METRICS ERROR:POST_DINUC\
            --ERROR_METRICS ERROR:CYCLE \
            --ERROR_METRICS ERROR:INSERT_LENGTH \
            --ERROR_METRICS ERROR:BASE_QUALITY \
            --ERROR_METRICS ERROR:MAPPING_QUALITY \
            --ERROR_METRICS ERROR:ONE_BASE_PADDED_CONTEXT
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        preemptible: preemptible
        disks: "local-disk " + disk_size + " HDD"
        cpu: 4
        memory: mem_gb + " GiB"
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
            }
    }
}

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    File wgs_coverage_interval_list
    Boolean use_fast_algorithm
    File reference_fasta
    File reference_index
    File reference_dict
    Int read_length

    String gatk_tag
    Int preemptible = 1
    Int mem_gb = 4
  }

  Float ref_size = size(reference_fasta, "GiB") + size(reference_fasta, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  String output_name = sub(basename(input_bam), "\.(bam|cram)$", "") + ".wgs_metrics"

  command {
    gatk \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{output_name} \
      USE_FAST_ALGORITHM=~{use_fast_algorithm} \
      READ_LENGTH=~{read_length}
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
    preemptible: preemptible
    memory: mem_gb + " GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "~{output_name}"
  }
}