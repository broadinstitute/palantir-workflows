version 1.0

workflow UMIAwareDuplicateMarking {
  input {
    File aligned_bam # aligned bam sorted by the query (read) name. 
    String output_basename
    Boolean remove_duplicates = false
    Boolean use_umi = true
  }

  if (use_umi){
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

    call SortSam as SortSamQueryName {
      input:
        input_bam = GroupByUMIs.grouped_bam,
        output_bam_basename = output_basename + ".grouped.queryname_sorted",
        sort_order = "queryname"
    }
  }

  File input_bam = if use_umi then select_first([SortSamQueryName.output_bam]) else aligned_bam

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
  }
}

task MarkDuplicates {
  input {
    File bam
    String output_basename
    Boolean use_umi
    Boolean remove_duplicates
  }

  String output_bam_basename = output_basename + ".duplicate_marked"
  Int disk_size = ceil(3 * size(bam, "GB")) + 128

  # We add the TAG_DUPLICATE_SET_MEMBERS flag for debugging/analysis purposes.
  # The flag should be removed in production to save storage cost.
  command <<<
    java -Xms8192m -jar /usr/gitc/picard.jar MarkDuplicates \
    -I ~{bam} \
    -O ~{output_bam_basename}.bam \
    --METRICS_FILE ~{output_basename}_duplicate_metrics.txt \ 
    --TAG_DUPLICATE_SET_MEMBERS \
    ~{true='--READ_ONE_BARCODE_TAG BX' false='' use_umi} \
    ~{true="--REMOVE_DUPLICATES" false="" remove_duplicates}

    samtools view -c -F 0x100 ~{output_bam_basename}.bam > duplicate_marked_read_count.txt
  >>>

  output {
    File duplicate_marked_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{output_basename}_duplicate_metrics.txt"
    Int duplciate_marked_read_count = read_int("duplicate_marked_read_count.txt")
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
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