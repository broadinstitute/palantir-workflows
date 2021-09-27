version 1.0

workflow UMIAwareDuplicateMarking {
  input {
    File aligned_bam # aligned bam sorted by the query (read) name. 
    String output_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  }

  # First sort the aligned bam by coordinate, so we can group duplicate sets using UMIs in the next step.
  call SortSam as SortSamFirst {
    input:
      input_bam = aligned_bam,
      output_bam_basename = output_basename + ".STAR_aligned.coorinate_sorted",
      sort_order = "coordinate"
  }

  # Further divide each duplicate set (a set of reads with the same insert start and end coordinates)
  # into subsets that share the same UMIs i.e. differenciate PCR duplicates from biological duplicates.
  # (biological duplicates are independent DNA molecules that are sheared such that the inserts are indistinguishable.)
  # input: a coordinate sorted bam
  # output: a coordinate sorted bam with UMIs (what are the generated tags?) .
  
  call GroupByUMIs {
    input:
      bam = SortSamFirst.output_bam,
      bam_index = select_first([SortSamFirst.output_bam_index, "bam_index_not_found"]),
      output_bam_basename = output_basename + ".grouped_by_UMI"
  }

  # input:
  # output: 
  call SortSam as SortSamQueryName {
    input:
      input_bam = GroupByUMIs.grouped_bam,
      output_bam_basename = output_basename + ".grouped.queryname_sorted",
      sort_order = "queryname"
  }

  call MarkDuplicates {
    input:
      bam = SortSamQueryName.output_bam,
      output_basename = output_basename
  }

  call SortSam as SortSamSecond {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename + ".duplicate_marked.coordinate_sorted.bam",
      sort_order = "coordinate"
  }

  call CollectMultipleMetrics {
    input:
      input_bam=SortSamSecond.output_bam,
      input_bam_index=select_first([SortSamSecond.output_bam_index, "bam_index_not_found"]),
      output_bam_prefix=output_basename,
      ref_dict=ref_dict,
      ref_fasta=ref_fasta,
      ref_fasta_index=ref_fasta_index,
      preemptible_tries=0
  }

  output {
    File duplicate_marked_bam = SortSamSecond.output_bam
    File duplicate_marked_bam_index = select_first([SortSamSecond.output_bam_index, "bam_index_not_found"])
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

task MarkDuplicates {
  input {
    File bam
    String output_basename
  }

  String output_bam_basename = output_basename + ".dulicate_marked.bam"

  Int disk_size = ceil(3 * size(bam, "GB")) + 128
  command <<<
    gatk MarkDuplicates -I ~{bam} --READ_ONE_BARCODE_TAG BX -O ~{output_bam_basename}.bam --METRICS_FILE ~{output_basename}.duplicate.metrics --ASSUME_SORT_ORDER queryname
  >>>

  output {
    File duplicate_marked_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{output_basename}.duplicate.metrics"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
  }
}

task SortSam {
  input {
    File input_bam
    String output_bam_basename
    String sort_order # "queryname" or "coordinate"
  }

  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 4.0
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 256

  String extra_args = if sort_order == "coordinate" then "CREATE_INDEX=true" else ""


  command {
    java -Xms8192m -jar /usr/picard/picard.jar \
    SortSam \
    INPUT=~{input_bam} \
    OUTPUT=~{output_bam_basename}.bam \
    SORT_ORDER=~{sort_order} \
    CREATE_MD5_FILE=true \
    MAX_RECORDS_IN_RAM=300000 \
    ~{extra_args}

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "16 GB"
    preemptible: 0
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File? output_bam_index = "~{output_bam_basename}.bai" # Create index for coordinate sort
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

task GroupByUMIs {
  input {
    File bam
    File bam_index
    String output_bam_basename
  }

  Int disk_space = ceil(2.2 * size(bam, "GB")) + 300
  command <<<
    umi_tools group -I ~{bam} --paired --no-sort-output --output-bam --stdout ~{output_bam_basename}.bam --umi-tag-delimiter "-" \
    --extract-umi-method tag --umi-tag RX --unmapped-reads use
  >>>

  output {
    File grouped_bam = "~{output_bam_basename}.bam"
  }

  runtime {
    docker : "us.gcr.io/tag-team-160914/tag-gtex-umi-tools:v1"
    disks : "local-disk " + disk_space + " HDD"
    preemptible: 0
    cpu: "8"
    memory: "64 GB" # Sato: is this too much?
  }
}

task CollectMultipleMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  File ref_flat = "gs://gcp-public-data--broad-references/hg38/v0/GRCh38_gencode.v27.refFlat.txt"

  command {
    java -Xms5000m -jar /usr/gitc/picard.jar CollectMultipleMetrics \
    INPUT=~{input_bam} \
    OUTPUT=~{output_bam_prefix} \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=CollectAlignmentSummaryMetrics \
    REFERENCE_SEQUENCE=~{ref_fasta}

    ls > ls.txt
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    memory: "8 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File ls = "ls.txt" # Take this out before release
    File alignment_summary_metrics = output_bam_prefix + ".alignment_summary_metrics"
    File insert_size_metrics = output_bam_prefix + ".insert_size_metrics"
    File insert_size_historgram = output_bam_prefix + ".insert_size_histogram.pdf"
    File base_distribution_by_cycle_metrics = output_bam_prefix + ".base_distribution_by_cycle_metrics"
    File quality_distribution_metrics = output_bam_prefix + ".quality_distribution_metrics"
  }
}