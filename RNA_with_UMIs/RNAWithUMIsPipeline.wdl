version 1.0

workflow RNAWithUMIsPipeline {
	input {
		File bam
		String read1Structure
		String read2Structure
		File starIndex
		String output_basename
		File gtf
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

	call SortSam as SortSamSTAR{
		input:
			input_bam = STAR.aligned_bam,
			output_bam_basename = "STAR.aligned.sorted"
	}

	call GroupByUMIs {
		input:
			bam = SortSamSTAR.output_bam,
			bam_index = SortSamSTAR.output_bam_index
	}

	call SortSamQuery {
		input:
			input_bam = GroupByUMIs.grouped_bam,
			output_bam_basename = "Grouped.queryname.sorted"
	}

	call MarkDuplicates {
		input:
			bam = SortSamQuery.output_bam
	}

	call SortSam {
		input:
			input_bam = MarkDuplicates.duplicate_marked_bam,
			output_bam_basename = output_basename
	}

	call GetSampleName {
		input:
			bam = bam
	}

	call RemoveUnwantedReads {
		input:
			bam = SortSam.output_bam
	}

	call rnaseqc2 {
		input:
			bam_file = SortSam.output_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name
	}

	call rnaseqc2 as rnaseqcRemoveDuplicateReads {
		input:
			bam_file = RemoveUnwantedReads.output_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name
	}

  output {
	File output_bam = SortSam.output_bam
	File output_bam_index = SortSam.output_bam_index
	File output_bam_md5 = SortSam.output_bam_md5
	File gene_tpm = rnaseqc2.gene_tpm
	File gene_counts = rnaseqc2.gene_counts
	File exon_counts = rnaseqc2.exon_counts
	File metrics = rnaseqc2.metrics
	File insertsize_distr = rnaseqc2.insertsize_distr
	File gene_dup = rnaseqc2.gene_dup
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
			--extract-umi-method tag --umi-tag RX
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

task STAR {
	input {
		File bam
		File starIndex
	}
	Int disk_space = ceil(2.2 * size(bam, "GB") + size(starIndex, "GB")) + 250

	command <<<
		echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
		mkdir star_index
		tar -xvvf ~{starIndex} -C star_index --strip-components=1

		STAR --readFilesIn ~{bam} --readFilesType SAM PE --readFilesCommand samtools view -h \
			--runMode alignReads --genomeDir star_index --outSAMtype BAM Unsorted --runThreadN 8 \
			--limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outSAMunmapped Within \
			--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 \
			--outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
			--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 --alignSoftClipAtReferenceEnds Yes --chimSegmentMin 15 --chimMainSegmentMultNmax 1 \
			--chimOutType WithinBAM SoftClip --chimOutJunctionFormat 0 --twopassMode Basic

	>>>

	runtime {
		docker : "us.gcr.io/tag-team-160914/neovax-tag-rnaseq:v1"
		disks : "local-disk " + disk_space + " HDD"
		memory : "64GB"
		cpu : "8"
		preemptible: 0
	}

	output {
		File aligned_bam = "Aligned.out.bam "
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
		gatk MarkDuplicates -I ~{bam} --READ_ONE_BARCODE_TAG BX -O duplicate.marked.bam --METRICS_FILE duplicate.metrics --ASSUME_SORT_ORDER queryname
	>>>

	output {
		File duplicate_marked_bam = "duplicate.marked.bam"
		File duplicate_metrics = "duplicate.metrics"
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

task rnaseqc2 {
	input {
		File bam_file
		File genes_gtf
		String sample_id
	}
	Int disk_space = ceil(size(bam_file, 'GB') + size(genes_gtf, 'GB')) + 100

	command {
		set -euo pipefail
		echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
		touch ~{sample_id}.fragmentSizes.txt
		touch ~{sample_id}.gene_duplicates.gct
		rnaseqc ~{genes_gtf} ~{bam_file} . -s ~{sample_id} -vv
		echo "  * compressing outputs"
		gzip *.gct
		echo $(date +"[%b %d %H:%M:%S] done")
	}

	output {
		File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
		File gene_counts = "${sample_id}.gene_reads.gct.gz"
		File exon_counts = "${sample_id}.exon_reads.gct.gz"
		File metrics = "${sample_id}.metrics.tsv"
		File insertsize_distr = "${sample_id}.fragmentSizes.txt"
		File gene_dup = "${sample_id}.gene_duplicates.gct.gz"
	}

	runtime {
		docker: "us.gcr.io/tag-team-160914/tag-gtex-rnaseqc-dup:v1"
		memory: "10GB"
		disks: "local-disk " + disk_space + " HDD"
		preemptible: 0
	}
}

task GetSampleName {
	input {
		File bam
	}

	parameter_meta {
		bam : {
			localization_optional : true
		}
	}

	command <<<
		gatk GetSampleName -I ~{bam} -O sample_name.txt
	>>>

	runtime {
		docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
		disks: "local-disk 100 HDD"
	}

	output {
		String sample_name = read_string("sample_name.txt")
	}
}

task RemoveUnwantedReads {
	input {
		File bam
	}

	parameter_meta {
		bam : {
			localization_optional : true
		}
	}

	command <<<
		gatk PrintReads --read-filter NotDuplicateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \
			--read-filter NotSecondaryAlignmentReadFilter --disable-tool-default-read-filters -I ~{bam} -O filtered.bam
	>>>

	runtime {
		docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
		disks: "local-disk 100 HDD"
	}

	output {
		File output_bam = "filtered.bam"
	}
}
