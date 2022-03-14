version 1.0

import "UMIAwareDuplicateMarking.wdl" as UmiMD

workflow RNAWithUMIsPipeline {
	input {
		File bam
		String read1Structure
		String read2Structure
		File starIndex
		String output_basename
		File gtf

		File ref
		File refIndex
		File refDict
		File refFlat
		File ribosomalIntervals

		File rnaseqc2_exon_bed
		Boolean use_umi
	}

	if (use_umi){
		call ExtractUMIs {
			input:
				bam = bam,
				read1Structure = read1Structure,
				read2Structure = read2Structure,
				output_prefix = output_basename
		}

		# We will not clip the non_use umi option, because that's really what we are comparing against.
		call SamToFastq {
			input:
				bam = ExtractUMIs.bam_umis_extracted,
				output_prefix = output_basename
		}

		# Adapter clipping
		call Fastp {
			input: 
				fastq1 = SamToFastq.fastq1,
				fastq2 = SamToFastq.fastq2,
				output_prefix = output_basename
		}

		call FastqToSam {
			input:
				fastq1=Fastp.fastq1_clipped,
				fastq2=Fastp.fastq2_clipped,
				sample_name=output_basename + "_clipped"
		}

		call FastQC as FastQCWithClipping {
			input:
				unmapped_bam = FastqToSam.unmapped_bam,
				sample_id = output_basename
		}
	}

	call FastQC {
		input:
			unmapped_bam = bam,
			sample_id = output_basename
	}

	File star_input_bam = if use_umi then select_first([FastqToSam.unmapped_bam]) else bam

	call STAR {
		input:
			bam = star_input_bam,
			starIndex = starIndex,
			transcriptome_ban = "IndelSoftclipSingleend"
	}

	# There is a STAR argument to do this, but this also gets the job done.
	call CopyReadGroupsToHeader {
		input:
			bam_with_readgroups = STAR.aligned_bam,
			bam_without_readgroups = STAR.transcriptome_bam
	}

	call UmiMD.UMIAwareDuplicateMarking {
		input:
			aligned_bam = STAR.aligned_bam,
			output_basename = output_basename,
			use_umi = use_umi,
			ubam = ExtractUMIs.bam_umis_extracted,
			ref_fasta = ref,
			ref_fasta_index = refIndex,
			ref_dict = refDict
	}

	call UmiMD.UMIAwareDuplicateMarking as UMIAwareDuplicateMarkingTranscriptome {
		input:
			aligned_bam = CopyReadGroupsToHeader.output_bam,
			output_basename = output_basename + "_transcriptome",
			use_umi = use_umi,
			remove_duplicates = true,
			ubam = ExtractUMIs.bam_umis_extracted,
			ref_fasta = ref,
			ref_fasta_index = refIndex,
			ref_dict = refDict
	}
	
	# This step should happen right after MD, before coordinate sorting.
	call FormatTranscriptomeUMI {
		input:
			prefix = output_basename + "_transcriptome_RSEM_formatted",
			input_bam = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_query_sorted_bam
	}


	# Extract the unaligned (i.e. unmapped) reads from the genome-aligned, duplicated marked bam,
	# which we will run FastQC on.
	call CreateUnalignedBam {
		input:
			input_bam = UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix = output_basename,
			preemptible_tries = 0
	}

	call FastQC as FastQCUnalignablereads {
		input:
			unmapped_bam = CreateUnalignedBam.unaligned_bam,
			sample_id = output_basename
	}

	# PLACEHOLDER FOR CROSSCHECK

	call GetSampleName {
		input:
			bam = bam
	}

	call rnaseqc2 {
		input:
			bam_file = UMIAwareDuplicateMarking.duplicate_marked_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name,
			exon_bed = rnaseqc2_exon_bed
	}

	call CollectRNASeqMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			ref_flat=refFlat,
			ribosomal_intervals=ribosomalIntervals,
			preemptible_tries=0
	}

	call CollectMultipleMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			preemptible_tries=0
	}

	call CollectMultipleMetrics as CollectMultipleMetricsWithoutDuplicateMarking {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			preemptible_tries=0
	}

	call CollectInsertSizeMetrics as InsertSizeTranscriptome {
		input:
			input_bam = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam,
			input_bam_index = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam_index,
			output_bam_prefix = GetSampleName.sample_name + "_transcriptome",
			preemptible_tries = 0
	}	

  output {
	File transcriptome_bam = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam
	File transcriptome_bam_index = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam_index
	File transcriptome_duplicate_metrics = UMIAwareDuplicateMarkingTranscriptome.duplicate_metrics
	File output_bam = UMIAwareDuplicateMarking.duplicate_marked_bam
	File output_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index
	File duplicate_metrics = UMIAwareDuplicateMarking.duplicate_metrics
	File gene_tpm = rnaseqc2.gene_tpm
	File gene_counts = rnaseqc2.gene_counts
	File exon_counts = rnaseqc2.exon_counts
	File rnaseqc2_metrics = rnaseqc2.metrics
	File rnaseqc2_insert_size_histogram = rnaseqc2.insert_size_histogram
	File fastqc_report = FastQC.fastqc_html
	File fastqc_table = FastQC.fastqc_data
	Float fastqc_adapter_content = FastQC.adapter_content
	File? fastqc_clipped_report = FastQCWithClipping.fastqc_html
	File? fastqc_clipped_table = FastQCWithClipping.fastqc_data
	Float? fastqc_clipped_adapter_content = FastQCWithClipping.adapter_content
	File? fastp_html = Fastp.html

	File genome_insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
	File transcriptome_insert_size_metrics = InsertSizeTranscriptome.insert_size_metrics
	File alignment_metrics = CollectMultipleMetrics.alignment_summary_metrics
	File rna_metrics = CollectRNASeqMetrics.rna_metrics
	Int pre_alignment_read_count = STAR.pre_alignment_read_count
	Int aligned_read_count = STAR.aligned_read_count
	Int transcriptome_read_count = STAR.transcriptome_read_count
	Float pct_reads_unmapped_mismatches = STAR.pct_reads_unmapped_mismatches
	Float pct_uniquely_mapped = STAR.pct_uniquely_mapped
	File formatted_transcriptome_bam = FormatTranscriptomeUMI.output_bam
	Int post_formatting_read_count = FormatTranscriptomeUMI.post_formatting_read_count
	
  }
}

task STAR {
	input {
		File bam
		File starIndex
		Int num_protrude_bases = 20
		String transcriptome_ban
	}

	Int disk_space = ceil(2.2 * size(bam, "GB") + size(starIndex, "GB")) + 250

	command <<<
		echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
		mkdir star_index
		tar -xvvf ~{starIndex} -C star_index --strip-components=1

		samtools view -c ~{bam} > pre_alignment_read_count.txt

		STAR --readFilesIn ~{bam} --readFilesType SAM PE --readFilesCommand samtools view -h \
			--runMode alignReads --genomeDir star_index --outSAMtype BAM Unsorted --runThreadN 8 \
			--outSAMunmapped Within \
			--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.2 \
			--outFilterMatchNminOverLread 0.2 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
			--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 --chimSegmentMin 15 --chimMainSegmentMultNmax 1 \
			--chimOutType WithinBAM SoftClip --chimOutJunctionFormat 0 --twopassMode Basic --quantMode TranscriptomeSAM \
			--quantTranscriptomeBan ~{transcriptome_ban} \
			--alignEndsProtrude ~{num_protrude_bases} ConcordantPair \

		samtools view -c -F 0x100 Aligned.out.bam > aligned_read_count.txt
		samtools view -c -F 0x100 Aligned.toTranscriptome.out.bam > transcriptome_read_count.txt
		ls > "ls.txt"

		grep "% of reads unmapped: too many mismatches" Log.final.out | cut -d "|" -f 2 | tr -d "[:space:]%" > pct_reads_unmapped_mismatches.txt
		grep "Uniquely mapped reads %" Log.final.out | cut -d "|" -f 2 | tr -d "[:space:]%" > pct_uniquely_mapped.txt
	>>>

	runtime {
		docker : "us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627"
		disks : "local-disk " + disk_space + " HDD"
		memory : "64GB"
		cpu : "8"
		preemptible: 0
	}

	output {
		File ls = "ls.txt"
		File aligned_bam = "Aligned.out.bam"
		File transcriptome_bam = "Aligned.toTranscriptome.out.bam"
		File splice_junction_table = "SJ.out.tab"
		File log_out = "Log.out"
		File star_metrics_log = "Log.final.out"
		Int pre_alignment_read_count = read_int("pre_alignment_read_count.txt")
		Int aligned_read_count = read_int("aligned_read_count.txt")
		Int transcriptome_read_count = read_int("transcriptome_read_count.txt")

		# STAR metrics
		Float pct_reads_unmapped_mismatches = read_float("pct_reads_unmapped_mismatches.txt")
		Float pct_uniquely_mapped = read_float("pct_uniquely_mapped.txt")
	}
}

task ExtractUMIs {
	input {
		File bam
		String read1Structure
		String read2Structure
		String output_prefix
	}

	Int disk_space = ceil(2.2 * size(bam, "GB")) + 50

	command <<<
		fgbio ExtractUmisFromBam --input ~{bam} \
			--read-structure ~{read1Structure} \
			--read-structure ~{read2Structure} \
			--molecular-index-tags RX \
			--output ~{output_prefix}_UMI_extracted.bam
	>>>

	runtime {
		docker : "quay.io/biocontainers/fgbio@sha256:a8e5cf58c318bffba3b2b694a3640ecd9e8106cee2e33b75710c0e8215138b6e"
		disks : "local-disk " + disk_space + " HDD"
		preemptible: 0
	}

	output {
		File bam_umis_extracted = "~{output_prefix}_UMI_extracted.bam"
	}
}

task rnaseqc2 {
	input {
		File bam_file
		File genes_gtf
		String sample_id
		File exon_bed
	}
	
	Int disk_space = ceil(size(bam_file, 'GB') + size(genes_gtf, 'GB')) + 100

	command {
		set -euo pipefail
		echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
		rnaseqc ~{genes_gtf} ~{bam_file} . -s ~{sample_id} -v --bed ~{exon_bed}
		echo "  * compressing outputs"
		gzip *.gct
		echo $(date +"[%b %d %H:%M:%S] done")
	}

	output {
		File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
		File gene_counts = "${sample_id}.gene_reads.gct.gz"
		File exon_counts = "${sample_id}.exon_reads.gct.gz"
		File insert_size_histogram = "${sample_id}.fragmentSizes.txt"
		File metrics = "${sample_id}.metrics.tsv"
	}

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/ckachulis/rnaseqc:2.4.2"
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

# STAR's transcriptome output does not include readgroups,
# so we simply append the read group lines from the header of the genome-aligned bam.
task CopyReadGroupsToHeader {
	input {
		File bam_with_readgroups
		File bam_without_readgroups
	}

	Int disk_size = ceil(2.0 * size([bam_with_readgroups, bam_without_readgroups], "GB")) + 10
	String basename = basename(bam_without_readgroups)

	command <<<
		samtools view -H ~{bam_without_readgroups} > header.sam
		samtools view -H ~{bam_with_readgroups} | grep "@RG" >> header.sam
		samtools reheader header.sam ~{bam_without_readgroups} > ~{basename}
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/samtools@sha256:0e49b0a5d91c203b8c07f5242277c2060b4b8ea54df8b1d123f990a1ad0588b2"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		File output_bam = basename
	}
}

task CollectRNASeqMetrics {
	input {
		File input_bam
		File input_bam_index
		String output_bam_prefix
		File ref_dict
		File ref_fasta
		File ref_fasta_index
		Int preemptible_tries
		File ref_flat
		File ribosomal_intervals
	}

	Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
	Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

	# This jar skips the header check of the ribosomal interval
	File picard_jar = "gs://broad-dsde-methods-takuto/RNA/picard_ribosomal_inserts.jar"
	File collapsed_ref_flat = "gs://broad-dsde-methods-takuto/RNA/resources/gencode_v34_UCSC_2_collapsed_refFlat.txt"

	command {
		java -Xms5000m -jar ~{picard_jar} CollectRnaSeqMetrics \
		REF_FLAT=~{ref_flat} \
		RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
		STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix}_rna_metrics.txt \
		RRNA_INS=~{output_bam_prefix}_rRNA_inserts.txt


		java -Xms5000m -jar ~{picard_jar} CollectRnaSeqMetrics \
		REF_FLAT=~{collapsed_ref_flat} \
		RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
		STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix}_rna_metrics_collapsed.txt \
		RRNA_INS=~{output_bam_prefix}_rRNA_inserts_collapsed.txt

		ls > ls.txt
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output {
		File rna_metrics = output_bam_prefix + "_rna_metrics.txt"
		File rna_metrics_collapsed = output_bam_prefix + "_rna_metrics_collapsed.txt"
		File rRNA_insert_size_histogram = output_bam_prefix + "_rRNA_inserts.txt"
		File ls = "ls.txt"
	}
}

#### Copied from warp/tasks/broad/Qc.wdl ####
#### Should be imported instead ####

task CrossCheckFingerprints {
	input {
		Array[File] input_bams
		Array[File] input_bam_indexes
		File haplotype_database_file
		String metrics_filename
		Float total_input_size
		Int preemptible_tries
		Float lod_threshold
		String cross_check_by
	}

	Int disk_size = ceil(total_input_size) + 20

	command {
		java -Dsamjdk.buffer_size=131072 \
		-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m \
		-jar /usr/picard/picard.jar \
		CrosscheckFingerprints \
		OUTPUT=~{metrics_filename} \
		HAPLOTYPE_MAP=~{haplotype_database_file} \
		EXPECT_ALL_GROUPS_TO_MATCH=true \
		INPUT=~{sep=' INPUT=' input_bams} \
		LOD_THRESHOLD=~{lod_threshold} \
		CROSSCHECK_BY=~{cross_check_by}
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
		preemptible: preemptible_tries
		memory: "3.5 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		File cross_check_fingerprints_metrics = "~{metrics_filename}"
	}
}


task FastQC {
	input {
		File unmapped_bam
		String sample_id
		Float? mem = 4
	}

	Int disk_size = ceil(size(unmapped_bam, "GiB") * 3)  + 100
	String bam_basename = basename(unmapped_bam, ".bam")

	command {
		perl /usr/tag/scripts/FastQC/fastqc ~{unmapped_bam} --extract -o ./
		mv ~{bam_basename}_fastqc/fastqc_data.txt ~{bam_basename}_fastqc_data.txt
		ls > ls.txt

		tail -n 2 ~{bam_basename}_fastqc_data.txt | head -n 1 | cut -f 2 > ~{bam_basename}_adapter_content.txt
	}
	
	runtime {
		docker: "us.gcr.io/tag-public/tag-tools:1.0.0"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem + "GB"
		cpu: "1"
	}

	output {
		File ls = "ls.txt"
		File fastqc_data = "~{bam_basename}_fastqc_data.txt"
		File fastqc_html = "~{bam_basename}_fastqc.html"
		Float adapter_content = read_float("~{bam_basename}_adapter_content.txt")
	}
}


# Unaligned reads after reads that failed to align to the reference,
# in contrast to reads that have yet to go through alignment.
task CreateUnalignedBam {
	input {
		File input_bam
		File input_bam_index
		String output_bam_prefix
		Int preemptible_tries
	}

	Int disk_size = ceil(size(input_bam, "GiB")) + 200

	command {
		samtools view -h -b -o ~{output_bam_prefix}_unaligned.bam ~{input_bam} '*'
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}

	output {
		File unaligned_bam = output_bam_prefix + "_unaligned.bam"
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

	command {
		java -Xms5000m -jar /usr/gitc/picard.jar CollectMultipleMetrics \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix} \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=CollectAlignmentSummaryMetrics \
		REFERENCE_SEQUENCE=~{ref_fasta}

		ls > ls.txt
		cp ~{output_bam_prefix}".alignment_summary_metrics" ~{output_bam_prefix}_alignment_summary_metrics.txt
		cp ~{output_bam_prefix}".insert_size_metrics" ~{output_bam_prefix}_insert_size_metrics.txt
	}	

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}

	output {
		File ls = "ls.txt"
		File alignment_summary_metrics = output_bam_prefix + "_alignment_summary_metrics.txt"
		File insert_size_metrics = output_bam_prefix + "_insert_size_metrics.txt"
	}
}

task CollectInsertSizeMetrics {
	input {
		File input_bam
		File input_bam_index
		String output_bam_prefix
		Int preemptible_tries
		File? interval_list
	}

	Int disk_size = ceil(size(input_bam, "GiB")) + 256

	command {
		java -Xms5000m -jar /usr/gitc/picard.jar CollectInsertSizeMetrics \
		I=~{input_bam} \
		O=~{output_bam_prefix}_insert_size_metrics.txt \
		H=~{output_bam_prefix}_insert_size_histogram.pdf \
		M=0.0
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "8 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output {
		File insert_size_metrics = output_bam_prefix + "_insert_size_metrics.txt"
		File insert_size_histogram = output_bam_prefix + "_insert_size_histogram.pdf"
	}
}

task SamToFastq {
	input {
		File bam
		String output_prefix
	}

	Int disk_size = 2*ceil(size(bam, "GiB")) + 128

	command {
		java -jar /usr/picard/picard.jar SamToFastq \
		I=~{bam} \
		FASTQ=~{output_prefix}_1.fastq \
		SECOND_END_FASTQ=~{output_prefix}_2.fastq

	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
		preemptible: 0
		memory: "3.5 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		File fastq1 = output_prefix + "_1.fastq"
		File fastq2 = output_prefix + "_2.fastq"
	}
}


task FastqToSam {
	input {
		File fastq1
		File fastq2
		String sample_name
	}

	Int disk_size = 3*ceil(size(fastq1, "GiB") + size(fastq2, "GiB")) + 256


	command {
		java -Xms8192m -jar /usr/picard/picard.jar FastqToSam \
		F1=~{fastq1} \
		F2=~{fastq2} \
		O=~{sample_name}.u.bam \
		SM=~{sample_name} \
		RG="RG1" \
		LB=~{sample_name} \
		PU="barcode1" \
		PL="ILLUMINA"
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
		disks: "local-disk " + disk_size + " HDD"
		cpu: "1"
		memory: "16 GB"
		preemptible: 0
	}

	output {
		File unmapped_bam = "~{sample_name}.u.bam"
	}
}

# Query sorts and...does other things.
task FormatTranscriptomeUMI {
  input {
    String prefix
    File input_bam
  }

  Int disk_gb = ceil(3*size(input_bam,"GB"))

  command {
    umi_tools prepare-for-rsem --tags UG,BX,RX \
      -I ~{input_bam} --stdout ~{prefix}.bam

    samtools view -c -F 0x100 ~{prefix}.bam > post_formatting_read_count.txt
  }
  
  output {
    File output_bam = "~{prefix}.bam"
    Int post_formatting_read_count = read_int("post_formatting_read_count.txt")
  }
  
  runtime {
    docker        : "us.gcr.io/tag-public/neovax-tag-rnaseq:v1"
    preemptible   : 0
    cpu           : "8"
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : "16GB"
  }
}

# for adapter clipping
task Fastp {
	input {
		File fastq1
		File fastq2
		String output_prefix
	}

	File adapter_fasta = "gs://broad-dsde-methods-takuto/RNA/resources/Illumina_adapters.fasta"
	Int disk_size = 5*ceil(size(fastq1, "GiB")) + 128

	command {
		fastp --in1 ~{fastq1} --in2 ~{fastq2} --out1 ~{output_prefix}_read1_trimmed.fastq --out2 ~{output_prefix}_read2_trimmed.fastq \
		--disable_quality_filtering --adapter_fasta ~{adapter_fasta} 
		ls > "ls.txt"
	}
	

	runtime {
		docker: "biocontainers/fastp:v0.20.1_cv1"
		memory: "8 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 0
	}

	output {
		File ls = "ls.txt"
		File fastq1_clipped = output_prefix + "_read1_trimmed.fastq"
		File fastq2_clipped = output_prefix + "_read2_trimmed.fastq"
		File html = "fastp.html"
		File json = "fastp.json"
	}

}

task AddNsToClippedReads {
	input {
		File fastq1
		File fastq2
		String output_prefix
		Int read_length = 151
	}

	Int disk_size = 5*ceil(size(fastq1, "GiB")) + 128
	File script = "gs://broad-dsde-methods-takuto/RNA/fix_read_length.py"

	command {
		python3 ~{script} --read_length ~{read_length} ~{fastq1} ~{output_prefix}_read1_padded.fastq
		python3 ~{script} --read_length ~{read_length} ~{fastq2} ~{output_prefix}_read2_padded.fastq
		ls > "ls.txt"
	}
	

	runtime {
		docker: "pegi3s/biopython"
		memory: "16 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 0
	}

	output {
		File ls = "ls.txt"
		File fastq1_padded = output_prefix + "_read1_padded.fastq"
		File fastq2_padded = output_prefix + "_read2_padded.fastq"
	}

}