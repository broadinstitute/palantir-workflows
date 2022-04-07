version 1.0

workflow UMIAwareDuplicateMarking {
  input {
    File aligned_bam # bam from aligner, not sorted.
    String output_basename
    Boolean remove_duplicates = false
    Boolean use_umi = true
    File? ubam # bam with UMIs extracted in RX tag
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  # Whether we use UMI or not, we need to sort by query name first.
  call SortSam as QueryNameSortAlignedBam {
    input:
      input_bam = aligned_bam,
      output_bam_basename = output_basename + ".queryname_sorted",
      sort_order = "queryname"
  }

  if (use_umi){
    # call MergeBamAlignment {
    #   input:
    #     aligned_bam = QueryNameSortAlignedBam.output_bam,
    #     ubam = select_first([ubam]), # if use_umi = true, then we should get the UMI-extracted ubam
    #     output_basename = output_basename,
    #     ref_fasta = ref_fasta,
    #     ref_fasta_index = ref_fasta_index,
    #     ref_dict = ref_dict
    # }

    # Recover RX tag
    call TransferReadTags {
      input:
        aligned_bam = QueryNameSortAlignedBam.output_bam,
        ubam = select_first([ubam]), # if use_umi = true, then we should get the UMI-extracted ubam
        output_basename = output_basename
    }

    # Sort the aligned bam by coordinate, so we can group duplicate sets using UMIs in the next step.
    call SortSam as SortSamFirst {
      input:
        input_bam = TransferReadTags.output_bam,
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

    call SortSam as SortSamQueryName {
      input:
        input_bam = GroupByUMIs.grouped_bam,
        output_bam_basename = output_basename + ".grouped.queryname_sorted",
        sort_order = "queryname"
    }
  }

  File input_bam = if use_umi then select_first([SortSamQueryName.output_bam]) else QueryNameSortAlignedBam.output_bam


  call MarkDuplicates {
    input:
      bam = input_bam,
      output_basename = output_basename,
      use_umi = use_umi,
      remove_duplicates = remove_duplicates
  }

  call SortSam as SortSamSecond {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename + ".duplicate_marked.coordinate_sorted",
      sort_order = "coordinate"
  }
  
  output {
    File duplicate_marked_query_sorted_bam = MarkDuplicates.duplicate_marked_bam
    File duplicate_marked_bam = SortSamSecond.output_bam
    File duplicate_marked_bam_index = select_first([SortSamSecond.output_bam_index, "bam_index_not_found"])
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    Int duplciate_marked_read_count = MarkDuplicates.duplciate_marked_read_count
    Int pre_transfer_count = TransferReadTags.pre_transfer_count
    Int post_transfer_count = TransferReadTags.post_transfer_count
  }
}

task TransferReadTags {
  input {
    File aligned_bam
    File ubam
    String output_basename
    File gatk_jar = "gs://broad-dsde-methods-takuto/RNA/gatk_transfer_read_tags.jar"
  }

  Int disk_size = ceil(2 * size(aligned_bam, "GB")) + ceil(2 * size(ubam, "GB")) + 128
  String output_bam_basename = output_basename + "_RX_transferred"
  
  command <<<
    java -jar ~{gatk_jar} TransferReadTags \
    -I ~{aligned_bam} \
    --unmapped-sam ~{ubam} \
    -O ~{output_bam_basename}.bam \
    --read-tags RX

    samtools view -c ~{aligned_bam} > pre_transfer_count.txt
    samtools view -c ~{output_bam_basename}.bam > post_transfer_count.txt
  >>>

  output {
    File output_bam = "~{output_bam_basename}.bam"
    Int pre_transfer_count = read_int("pre_transfer_count.txt")
    Int post_transfer_count = read_int("post_transfer_count.txt")
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
  }
}

task MergeBamAlignment {
  input {
    File picard_jar = "gs://broad-dsde-methods-takuto/RNA/picard_skip_header_check.jar"
    File aligned_bam
    File ubam
    String output_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  Int disk_size = ceil(2 * size(aligned_bam, "GB")) + ceil(2 * size(ubam, "GB")) + 128

  String output_bam_basename = output_basename + "_merged"
  command <<<
    java -Xms8192m -jar ~{picard_jar} MergeBamAlignment \
    UNMAPPED_BAM=~{ubam} \
    ALIGNED_BAM=~{aligned_bam} \
    OUTPUT=~{output_bam_basename}.bam \
    REFERENCE_SEQUENCE=~{ref_fasta} \
    SORT_ORDER=queryname \
    MATCHING_DICTIONARY_TAGS=null
  >>>

  output {
    File merged_bam = "~{output_bam_basename}.bam"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
  }
}

task MarkDuplicates {
  input {
    File bam
    String output_basename
    Boolean use_umi
    Boolean remove_duplicates
    File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
  }

  String output_bam_basename = output_basename + "_duplicate_marked"
  Int disk_size = ceil(3 * size(bam, "GB")) + 128

  # We add the TAG_DUPLICATE_SET_MEMBERS flag for debugging/analysis purposes.
  # The flag should be removed in production to save storage cost.
  command <<<
    bash ~{monitoring_script} > monitoring.log &

    java -Xms12g -jar /usr/gitc/picard.jar MarkDuplicates \
    INPUT=~{bam} \
    OUTPUT=~{output_bam_basename}.bam \
    METRICS_FILE=~{output_basename}_duplicate_metrics.txt \
    TAG_DUPLICATE_SET_MEMBERS=true \
    ~{true="READ_ONE_BARCODE_TAG=BX" false="" use_umi} \
    ~{true="REMOVE_DUPLICATES=true" false="" remove_duplicates}

    samtools view -c -F 0x100 ~{output_bam_basename}.bam > duplicate_marked_read_count.txt

    ls > ls.txt
    ls /usr >> ls.txt
  >>>

  output {
    File ls = "ls.txt"
    File duplicate_marked_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{output_basename}_duplicate_metrics.txt"
    Int duplciate_marked_read_count = read_int("duplicate_marked_read_count.txt")
    File monitoring_log = "monitoring.log"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z" # Docker needs to have MarkDuplicates and samtools. "samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438" gives permission denier for docker pull.
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
    preemptible: 0
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
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
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

  Int disk_space = ceil(4 * size(bam, "GB")) + 300
  command <<<
    umi_tools group -I ~{bam} --paired --no-sort-output --output-bam -S ~{output_bam_basename}.bam --umi-tag-delimiter "-" \
    --extract-umi-method tag --umi-tag RX --unmapped-reads use \
    --group-out ~{output_bam_basename}_umi_groups.txt
  >>>

  output {
    File grouped_bam = "~{output_bam_basename}.bam"
    File groups_file = "~{output_bam_basename}_umi_groups.txt"
  }

  runtime {
    docker : "us.gcr.io/tag-team-160914/tag-gtex-umi-tools:v1"
    disks : "local-disk " + disk_space + " HDD"
    preemptible: 0
    cpu: "8"
    memory: "64 GB" # Sato: is this too much?
  }
}