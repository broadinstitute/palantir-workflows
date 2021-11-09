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

		File contamination_sites_ud = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.UD"
		File contamination_sites_bed = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.bed"
		File contamination_sites_mu = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.contam.mu"
	}

	call ExtractUMIs {
		input:
			bam = bam,
			read1Structure = read1Structure,
			read2Structure = read2Structure
	}

	call FastQC {
		input:
			unmapped_bam = bam,
			sample_id = output_basename
	}

	call STAR {
		input:
			bam = ExtractUMIs.bam_umis_extracted,
			starIndex = starIndex
	}

	call CopyReadGroupsToHeader {
		input:
			bam_with_readgroups = STAR.aligned_bam,
			bam_without_readgroups = STAR.transcriptome_bam
	}

	call UmiMD.UMIAwareDuplicateMarking {
		input:
			aligned_bam = STAR.aligned_bam,
			output_basename = output_basename
	}

	call UmiMD.UMIAwareDuplicateMarking as UMIAwareDuplicateMarkingTranscriptome {
		input:
			aligned_bam = CopyReadGroupsToHeader.output_bam,
			output_basename = output_basename + ".transcriptome"
	}

	# TODO :Add crosscheck here.
	# call CrossCheckFingerprints {
	# 	input:
	# 		input_bams = [UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam],
	# 		input_bam_indexes = [UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam_index],
	# 		haplotype_database_file = ,
	# 		metrics_filename = 
	# 		total_input_size = 
	# 		preemptible_tries = 0,
	# 		lod_threshold = 0,
	# 		cross_check_by = 
	# }

	call CheckContamination {
		input:
			input_bam = UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			contamination_sites_ud = contamination_sites_ud,
			contamination_sites_bed = contamination_sites_bed,
			contamination_sites_mu = contamination_sites_mu,
			ref_fasta = ref,
			ref_fasta_index = refIndex,
			output_prefix = output_basename,
			preemptible_tries = 0
	}

	call GetSampleName {
		input:
			bam = bam
	}

	call rnaseqc2 {
		input:
			bam_file = UMIAwareDuplicateMarking.duplicate_marked_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name
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
	File metrics = rnaseqc2.metrics
		
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
			--chimOutType WithinBAM SoftClip --chimOutJunctionFormat 0 --twopassMode Basic --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend

	>>>

	runtime {
		docker : "us.gcr.io/tag-team-160914/neovax-tag-rnaseq:v1"
		disks : "local-disk " + disk_space + " HDD"
		memory : "64GB"
		cpu : "8"
		preemptible: 0
	}

	output {
		File aligned_bam = "Aligned.out.bam"
		File transcriptome_bam = "Aligned.toTranscriptome.out.bam"
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

task rnaseqc2 {
	input {
		File bam_file
		File genes_gtf
		String sample_id
	}
	Int disk_space = ceil(size(bam_file, 'GB') + size(genes_gtf, 'GB')) + 100
	File exon_bed = "gs://gtex-resources/GENCODE/gencode.v26.GRCh38.insert_size_intervals_geq1000bp.bed"

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
	File picard_jar = "gs://broad-dsde-methods-takuto/hydro.gen/picard_ignore_ribosomal_header.jar"

	command {
		java -Xms5000m -jar ~{picard_jar} CollectRnaSeqMetrics \
		REF_FLAT=~{ref_flat} \
		RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
		STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix}.rna_metrics
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output {
		File rna_metrics = output_bam_prefix + ".rna_metrics"
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
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}

	output {
		File ls = "ls.txt"
		File alignment_summary_metrics = output_bam_prefix + ".alignment_summary_metrics"
		File insert_size_metrics = output_bam_prefix + ".insert_size_metrics"
	}
}

task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Boolean disable_sanity_check = false
  }

  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"]))
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "7.5 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.0-c1cba76e979904eb69c31520a0d7f5be63c72253-1626442707"
    cpu: 2
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}

