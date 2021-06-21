version 1.0

workflow RNAWithUMIsPipeline {
	input {
		File bam
		String read1Structure
		String read2Structure
		File starIndex
		String output_basename
	}

	call ExtractUMIs {
		input:
			bam = bam,
			read1Structure = read1Structure,
			read2Structure = read2Structure
	}

	call STAR {
		input:
			bam = ExtractUMIs.bam_umis_extracted,
			starIndex = starIndex
	}

	call GroupByUMIs {
		input:
			bam = STAR.aligned_bam
	}

	call MarkDuplicates {
		input:
			bam = GroupByUMIs.grouped_bam
	}

	call SortSam {
		input:
			input_bam = MarkDuplicates.duplicate_marked_bam,
			output_bam_basename = output_basename
	}

  output {
	File output_bam = SortSam.output_bam
	File output_bam_index = SortSam.output_bam_index
	File output_bam_md5 = SortSam.output_bam_md5
  }
}

task GroupByUMIs {
	input {
		File bam
	}

	Int disk_space = ceil(2.2 * size(bam, "GB")) + 50
	command <<<
		umi_tools group -I ~{bam} --paired --no-sort-output --unpaired-reads output --output-bam --stdout umis.grouped.bam --umi-tag-delimiter "-"
	>>>

	output {
		File grouped_bam = "umis.grouped.bam"
	}

	runtime {
		docker : "us.gcr.io/tag-team-160914/tag-gtex-umi-tools:v1"
		disks : "local-disk " + disk_space + " HDD"
		preemptible: 0
	}
}

task STAR {
	input {
		File bam
		File starIndex
	}
	Int disk_space = ceil(2.2 * size(bam, "GB") + size(starIndex, "GB")) + 50

	command <<<
		echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
		mkdir star_index
		tar -xvvf ${star_index} -C star_index --strip-components=1

		STAR --inputBAMFile ~{bam} --runMode inputAlignmentsFromBam --genomeDir star_index --outSAMtype BAM > star.aligned.bam
	>>>

	runtime {
		docker : "us.gcr.io/tag-team-160914/neovax-tag-rnaseq:v1"
		disks : "local-disk " + disk_space + " HDD"
		preemptible: 0
	}

	output {
		File aligned_bam = "star.aligned.bam"
	}
}

task ExtractUMIs {
	input {
		File bam
		String read1Structure
		String read2Structure
	}

	Int disk_space = ceil(2.2 * size(bam, "GB")) + 50

	command <<<
		fgbio ExtractUmisFromBam --input ~{bam} \
			--read-structure ~{read1Structure} \
			--read-structure ~{read2Structure} \
			--molecular-index-tags RX \
			--output extractUMIs.out.bam
	>>>

	runtime {
		docker : "quay.io/biocontainers/fgbio@sha256:a8e5cf58c318bffba3b2b694a3640ecd9e8106cee2e33b75710c0e8215138b6e"
		disks : "local-disk " + disk_space + " HDD"
		preemptible: 0
	}

	output {
		File bam_umis_extracted = "extractUMIs.out.bam"
	}
}

task MarkDuplicates {
	input {
		File bam
	}

	Int disk_size = ceil(2.2 * size(bam, "GB")) + 50
	command <<<
		gatk MarkDuplicates -I ~{bam} --BARCODE_TAG BX --DUPLEX_UMI -O duplicate.marked.bam
	>>>

	output {
		File duplicate_marked_bam = "duplicate.marked.bam"
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