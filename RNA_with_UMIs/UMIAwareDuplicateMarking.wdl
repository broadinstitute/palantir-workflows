version 1.0

# The UMI-*Un*aware version
workflow UMIAwareDuplicateMarking {
  input {
    File aligned_bam
    String output_basename
  }

#  call SortSam as SortSamFirst {
#    input:
#      input_bam = aligned_bam,
#      output_bam_basename = "STAR.aligned.sorted"
#  }

#  call GroupByUMIs {
#    input:
#      bam = SortSamFirst.output_bam,
#      bam_index = SortSamFirst.output_bam_index
#  }

#  call SortSamQuery {
#    input:
#      input_bam = GroupByUMIs.grouped_bam,
#      output_bam_basename = "Grouped.queryname.sorted"
#  }

  call MarkDuplicates {
    input:
      bam = aligned_bam
  }

  call SortSam {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename
  }

  output {
    File duplicate_marked_bam = SortSam.output_bam
    File duplicate_marked_bam_index = SortSam.output_bam_index
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

task MarkDuplicates {
  input {
    File bam
  }

  String basename = basename(bam, ".bam")

  Int disk_size = ceil(2.2 * size(bam, "GB")) + 50
  command <<<
    gatk MarkDuplicates -I ~{bam} --READ_ONE_BARCODE_TAG BX -O ~{basename}.duplicate.marked.bam --METRICS_FILE ~{basename}.duplicate.metrics --ASSUME_SORT_ORDER queryname
  >>>

  output {
    File duplicate_marked_bam = "~{basename}.duplicate.marked.bam"
    File duplicate_metrics = "~{basename}.duplicate.metrics"
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
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 3.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

  command {
    java -Xms4000m -jar /usr/picard/picard.jar \
    SortSam \
    INPUT=~{input_bam} \
    OUTPUT=~{output_bam_basename}.bam \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true \
    MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MiB"
    preemptible: 0
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

task SortSamQuery {
  input {
    File input_bam
    String output_bam_basename
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 3.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

  command {
    java -Xms4000m -jar /usr/picard/picard.jar \
    SortSam \
    INPUT=~{input_bam} \
    OUTPUT=~{output_bam_basename}.bam \
    SORT_ORDER="queryname" \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true \
    MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MiB"
    preemptible: 0
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task GroupByUMIs {
  input {
    File bam
    File bam_index
  }

  Int disk_space = ceil(2.2 * size(bam, "GB")) + 300
  command <<<
    umi_tools group -I ~{bam} --paired --no-sort-output --output-bam --stdout umis.grouped.bam --umi-tag-delimiter "-" \
    --extract-umi-method tag --umi-tag RX --unmapped-reads use
  >>>

  output {
    File grouped_bam = "umis.grouped.bam"
  }

  runtime {
    docker : "us.gcr.io/tag-team-160914/tag-gtex-umi-tools:v1"
    disks : "local-disk " + disk_space + " HDD"
    preemptible: 0
    cpu: "8"
    memory: "52GB"
  }
}