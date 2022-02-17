	### Clip adapter reads experiment ####
	### Disabled and moved here since it didn't improve the performance much ####
	# Peprocessing for clipping by fastp
	call SamToFastq {
		input:
			bam = bam,
			output_prefix = output_basename
	}

	# Adapter clipping
	call Fastp {
		input: 
			fastq1 = SamToFastq.fastq1,
			fastq2 = SamToFastq.fastq2,
			output_prefix = output_basename
	}

	call AddNsToClippedReads {
		input:
			fastq1 = Fastp.fastq1_clipped,
			fastq2 = Fastp.fastq2_clipped,
			output_prefix = output_basename
	}

	call FastqToSam {
		input:
			fastq1=AddNsToClippedReads.fastq1_padded,
			fastq2=AddNsToClippedReads.fastq2_padded,
			sample_name=output_basename + "_clipped_padded"
	}

	call ExtractUMIs as ExtractUMIsClipped {
		input:
			bam = FastqToSam.unmapped_bam,
			read1Structure = read1Structure,
			read2Structure = read2Structure
	}

	call STAR as STARClipped {
		input:
			bam = ExtractUMIsClipped.bam_umis_extracted,
			starIndex = starIndex
	}

	call UmiMD.UMIAwareDuplicateMarking as UMIAwareDuplicateMarkingClipped {
		input:
			aligned_bam = STARClipped.aligned_bam,
			output_basename = output_basename + "_clipped"
	}

	call CollectRNASeqMetrics as CollectRNASeqMetricsClipped {
		input:
			input_bam=UMIAwareDuplicateMarkingClipped.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarkingClipped.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name + "_clipped",
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			ref_flat=refFlat,
			ribosomal_intervals=ribosomalIntervals,
			preemptible_tries=0
	}

	call CollectMultipleMetrics as CollectMultipleMetricsClipped {
		input:
			input_bam=UMIAwareDuplicateMarkingClipped.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarkingClipped.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name + "_clipped",
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			preemptible_tries=0
	}


# for adapter clipping
task Fastp {
	input {
		File fastq1
		File fastq2
		String output_prefix
	}

	Int disk_size = 5*ceil(size(fastq1, "GiB")) + 128

	command {
		fastp --in1 ~{fastq1} --in2 ~{fastq2} --out1 ~{output_prefix}_read1_trimmed.fastq --out2 ~{output_prefix}_read2_trimmed.fastq
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
		memory: "8 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 0
	}

	output {
		File ls = "ls.txt"
		File fastq1_padded = output_prefix + "_read1_padded.fastq"
		File fastq2_padded = output_prefix + "_read2_padded.fastq"
	}

}