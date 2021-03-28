version 1.0

import "../Imputation.wdl" as Imputation
import "../Structs.wdl" as structs

workflow validateImputation {
	input {
		File validationArrays
		File validationArraysIndex
		File validationWGS

		File? af_resource
		String? gnomad_vcf_prefix
		String? chr_prefix
		Array[String]? expressions

		File? reference_fasta
		File? reference_index
		File? reference_dict

		File subpopulation_af_expression
		File? sample_map

		Array[ReferencePanelContig] referencePanelContigs
	}

	call Imputation.ImputationPipeline {
		input:
			multi_sample_vcf = validationArrays,
			multi_sample_vcf_index = validationArraysIndex,
			referencePanelContigs = referencePanelContigs,
			perform_extra_qc_steps = false
	}

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

	call PearsonCorrelation {
		input:
			evalVcf = ImputationPipeline.imputed_multisample_vcf,
			truthVcf = validationWGS,
			af_expressions = subpopulation_af_expression,
			sample_map = sample_map,
			af_resource = select_first([af_resource, GatherVCFsCloud.vcf_out]),
			output_basename = "validation",
			missingIsHomRef = true
	}

	call plotCorrelations {
		input :
			correlations = PearsonCorrelation.correlations
	}

	output {
		File correlations = PearsonCorrelation.correlations
		File accuracy = PearsonCorrelation.accuracy
		File correlations_plot = plotCorrelations.correlations_plot
		File aggregated_imputation_metrics = ImputationPipeline.aggregated_imputation_metrics

		File imputed_multisample_vcf = ImputationPipeline.imputed_multisample_vcf
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
		Boolean missingIsHomRef = false
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
		-OA ~{output_basename}.accuracy.tsv ~{"-nBins " + n_bins} ~{"-firstBinRightEdge " + right_edge_first_bin} ~{"-minAfForAccuracyMetrics " + min_af_for_accuracy_metrics} \
		~{if missingIsHomRef then "-missingIsHomRef" else ""}
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:59b5115c3c0a521ff69adfb8b5cf7a826bea52155aef3918a3faa6914cd5d535"
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
	}

	command <<<
		Rscript -<< "EOF"
		library(readr)
		library(dplyr)
		library(ggplot2)

		corr <- read_tsv("~{correlations}")
		corr_gathered <- gather(correlations, key="type", value="correlation", snp_correlation, indel_correlation) %>%
				mutate(variant_type=ifelse(type=="snp_correlation", "snp", "indel"))
		ggplot(corr_gathered %>% filter(!is.na(correlations)), aes(aes(x=bin_center, y=correlation^2)) +
			geom_point(size=0.2, alpha=0.1) +
			geom_smooth(se=FALSE) +
			facet_grid(variant_type~af_annotation) + scale_x_log10() + theme_bw() +
			xlab("Minor Allele Frequency") + ylab(bquote(R^2))

		ggsave(filename="correlation_plot.png")

		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File correlations_plot = "correlation_plot.png"
	}
}