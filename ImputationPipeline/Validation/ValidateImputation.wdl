version 1.0

import "../Imputation.wdl" as Imputation
import "../Structs.wdl" as structs
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/ImputationPipeline/Imputation.wdl" as ImputationMain

workflow validateImputation {
	input {
		File validationArrays #input array vcf to use for validation
		File validationArraysIndex
		File validationWGS #input wgs vcf to use for validation
		String branch #name of branch being tested.  used for display in plots, does not effect computation

		File? af_resource #resource vcf containing population allele frequencies for sites
		String? gnomad_vcf_prefix #gnomad vcf prefix, for annotating with population af if af_resource not provided
		String? chr_prefix #"chr" is hg38, "" is hg19
		Array[String]? expressions #population allele frequency annotation to extract from gnomad if af_resource is not provided

		File? reference_fasta
		File? reference_index
		File? reference_dict

		File subpopulation_af_expression #File which maps sample to what allele frequence annotation should be used to extract subpopulation af from af_reource for each sample.  ":" used as separator
		File samplePopulations #tsv which maps sample to it's subpopulation, for grouping samples in plots
		File? sample_map #File which maps sample names in array to sample names in wgs.  ":" used as separator

		Array[ReferencePanelContig] referencePanelContigs
		File haplotype_database
	}

	#run imputation on this branch
	call Imputation.ImputationPipeline {
		input:
			multi_sample_vcf = validationArrays,
			multi_sample_vcf_index = validationArraysIndex,
			referencePanelContigs = referencePanelContigs,
			perform_extra_qc_steps = false,
			haplotype_database = haplotype_database
	}

	#run imputation on main branch
	call ImputationMain.ImputationPipeline as ImputationPipelineMain {
		input:
			multi_sample_vcf = validationArrays,
			multi_sample_vcf_index = validationArraysIndex,
			referencePanelContigs = referencePanelContigs,
			perform_extra_qc_steps = false,
			haplotype_database = haplotype_database
	}

	#if we do not have a af_resource from input, we need to make one by adding gnomad allele frequencies to a sites-only vcf of the imputed sites
	if(!defined(af_resource)) {
		call MakeSitesOnlyVcf {
			input:
				vcf = ImputationPipeline.imputed_multisample_vcf,
				vcfIndex = ImputationPipeline.imputed_multisample_vcf_index
		}

		scatter (chrom_i in range(22)) {
		  Int chrom = chrom_i + 1 ## range(22) gives you 0..21

		  call AnnotateWithAF_t {
			input:
			  vcf = MakeSitesOnlyVcf.sites_only_vcf,
			  vcf_index  = MakeSitesOnlyVcf.sites_only_vcf_index,
			  gnomad_vcf = select_first([gnomad_vcf_prefix]) + chrom + ".vcf.bgz",
			  gnomad_vcf_index = select_first([gnomad_vcf_prefix]) + chrom + ".vcf.bgz.tbi",
			  mem = 32, # need more memory for chrom 1
			  interval = select_first([chr_prefix]) + chrom,
			  ref_fasta = select_first([reference_fasta]),
			  ref_fasta_index = select_first([reference_index]),
			  ref_dict = select_first([reference_dict]),
			  output_basename = "annotated.snps_only." + chrom,
			  expressions = select_first([expressions])
		  }
		}

		call AnnotateWithAF_t as AnnotateX {
			input:
				vcf = MakeSitesOnlyVcf.sites_only_vcf,
				vcf_index  = MakeSitesOnlyVcf.sites_only_vcf,
				gnomad_vcf = select_first([gnomad_vcf_prefix]) + "X.vcf.bgz",
				gnomad_vcf_index = select_first([gnomad_vcf_prefix]) + "X.vcf.bgz.tbi",
				mem = 32,
				interval = select_first([chr_prefix]) + "X",
				ref_fasta = select_first([reference_fasta]),
				ref_fasta_index = select_first([reference_index]),
				ref_dict = select_first([reference_dict]),
				output_basename = "annotated.snps_only.X",
				expressions = select_first([expressions])
		}

		call GatherVCFsCloud {
			input:
				vcfs = flatten([AnnotateWithAF_t.annotated_vcf, [AnnotateX.annotated_vcf]])
		}
	}

	#imputation of this branch
	call PearsonCorrelation {
		input:
			evalVcf = ImputationPipeline.imputed_multisample_vcf,
			truthVcf = validationWGS,
			af_expressions = subpopulation_af_expression,
			sample_map = sample_map,
			af_resource = select_first([af_resource, GatherVCFsCloud.vcf_out]),
			output_basename = "validation"
	}

	#imputation performance of main branch
	call PearsonCorrelation as PearsonCorrelationMain{
		input:
			evalVcf = ImputationPipelineMain.imputed_multisample_vcf,
			truthVcf = validationWGS,
			af_expressions = subpopulation_af_expression,
			sample_map = sample_map,
			af_resource = select_first([af_resource, GatherVCFsCloud.vcf_out]),
			output_basename = "validation"
	}

	call plotCorrelations {
		input :
			correlations = PearsonCorrelation.correlations,
			correlationsMain = PearsonCorrelationMain.correlations,
			samplePopulations = samplePopulations,
			branch = branch
	}

	output {
		File correlations = PearsonCorrelation.correlations
		File accuracy = PearsonCorrelation.accuracy
		File correlations_plot = plotCorrelations.correlations_plot
		File aggregated_imputation_metrics = ImputationPipeline.aggregated_imputation_metrics

		File imputed_multisample_vcf = ImputationPipeline.imputed_multisample_vcf
		File imputed_multisample_vcf_main = ImputationPipelineMain.imputed_multisample_vcf
	}
}


task MakeSitesOnlyVcf {
	input {
		File vcf
		File vcfIndex
	}

	Int disk_size = 100 + ceil(2.2*size(vcf,"GB"))

	command <<<
		gatk MakeSitesOnlyVcf -I ~{vcf} -O sites_only.vcf.gz
	>>>

	runtime {
		docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}

	output {
		File sites_only_vcf = "sites_only.vcf.gz"
		File sites_only_vcf_index = "sites_only.vcf.gz.tbi"
	}
}

task AnnotateWithAF_t {
  input {
    File vcf
    File vcf_index
    File gnomad_vcf
    File gnomad_vcf_index
    String interval
    File ref_fasta
    File ref_fasta_index
    String output_basename
    File ref_dict
    Int mem = 16
    Array[String] expressions
  }

  Int disk_size = 400

  command <<<
    gatk VariantAnnotator -R ~{ref_fasta} -V ~{vcf} -L ~{interval} -O ~{output_basename}.vcf.gz  \
    --resource:gnomad ~{gnomad_vcf}  --expression ~{sep=" --expression " expressions} -LE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB" # some of the gnomad vcfs are like 38 gigs so maybe need more ?
  }
  output {
    File annotated_vcf = "~{output_basename}.vcf.gz"
    File annotated_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task GatherVCFsCloud {
	input {
		Array[File] vcfs
	}

	Int disk_size = ceil(2* size(vcfs, "GB")) + 10

	parameter_meta {
		vcfs : {
			localization_optional : true
		}
	}

	command <<<
		gatk GatherVcfsCloud -I ~{sep=" -I " vcfs} --gather-type CONVENTIONAL -O merged.vcf.gz
	>>>

	runtime {
    		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    		disks: "local-disk " + disk_size + " HDD"
    		memory: "16 GB"
    	}

	output {
		File vcf_out = "merged.vcf.gz"
		File vcf_out_index = "merged.vcf.gz.tbi"
	}
}

task PearsonCorrelation {
	input {
		File evalVcf
		File truthVcf
		File af_expressions
		File? sample_map
		String output_basename
		File af_resource
		File? sites
		String? intervals
		String? dosage_field
		Int? n_bins
		Float? right_edge_first_bin
		Float? min_af_for_accuracy_metrics
		Int mem = 16
	}

	parameter_meta {
		evalVcf : {
			localization_optional : true
		}
		truthVcf : {
			localization_optional : true
		}
		af_resource : {
			localization_optional : true
		}
	}

	command <<<
		set -xeuo pipefail


		cp ~{af_expressions} af_expressions.list
		touch sample_map.list
		~{"cp " + sample_map + " sample_map.list"}

		gatk --java-options "-Xmx~{mem - 2}G" ArrayImputationCorrelation --eval ~{evalVcf} --truth ~{truthVcf} --af-annotations af_expressions.list --resource ~{af_resource} \
		~{"--ids " + sites} ~{"-L " + intervals} --sample-map sample_map.list ~{"--dosage-field " + dosage_field} -O ~{output_basename}.correlations.tsv \
		-OA ~{output_basename}.accuracy.tsv ~{"-nBins " + n_bins} ~{"-firstBinRightEdge " + right_edge_first_bin} ~{"-minAfForAccuracyMetrics " + min_af_for_accuracy_metrics}
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:959899a92957d7987033e420a7b8b3469b7d9456023383725aa16215b6e75318"
		disks: "local-disk 100 HDD"
		memory: mem + " GB"
	}

	output {
		File correlations = "~{output_basename}.correlations.tsv"
		File accuracy = "~{output_basename}.accuracy.tsv"
	}
}

task plotCorrelations {
	input {
		File correlations
		File correlationsMain
		File samplePopulations
		String branch
	}

	command <<<
		Rscript -<< "EOF"
		library(readr)
		library(dplyr)
		library(ggplot2)
		library(tidyr)

		corr <- read_tsv("~{correlations}") %>% mutate(branch="~{branch}")
		corr_main <- read_tsv("~{correlationsMain}") %>% mutate(branch="main")
		populations <- read_tsv("~{samplePopulations}")
		corr <- bind_rows(corr, corr_main)
		corr <- inner_join(corr,populations)
		corr_gathered <- gather(corr, key="type", value="correlation", snp_correlation, indel_correlation) %>%
				mutate(variant_type=ifelse(type=="snp_correlation", "snp", "indel"))
		ggplot(corr_gathered %>% filter(!is.na(correlation)), aes(x=bin_center, y=correlation^2)) +
			geom_point(size=0.2, alpha=0.1, aes(color=branch)) +
			geom_smooth(se=FALSE, aes(color=branch)) +
			facet_grid(variant_type~sub_population) + scale_x_log10() + theme_bw() +
			xlab("Minor Allele Frequency") + ylab(bquote(R^2))

		ggsave(filename="correlation_plot.png")

		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse@sha256:f9671fa9329160cc57f76ec822896467f420f4b8d55bf3a811293b9f94283a3e"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File correlations_plot = "correlation_plot.png"
	}
}