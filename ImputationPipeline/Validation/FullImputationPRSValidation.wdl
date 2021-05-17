version 1.0

import "ValidateImputation.wdl" as ValidateImputation
import "ValidateScoring.wdl" as ValidateScoring

workflow FullImputationPRSValidation {
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

		String population_basename

		File population_loadings
		File population_meansd
		File population_pcs
		File pruning_sites_for_pca # and the sites used for PCA
		File population_vcf

		File weights
		File sample_name_map

		String branch

		Int wgs_vcf_to_plink_mem = 8
		File haplotype_database
	}

	call ValidateImputation.validateImputation {
		input:
			validationArrays = validationArrays,
			validationArraysIndex = validationArraysIndex,
			validationWGS = validationWGS,
			af_resource = af_resource,
			gnomad_vcf_prefix = gnomad_vcf_prefix,
			chr_prefix = chr_prefix,
			expressions = expressions,
			reference_fasta = reference_fasta,
			reference_index = reference_index,
			reference_dict = reference_dict,
			subpopulation_af_expression = subpopulation_af_expression,
			sample_map = sample_map,
			referencePanelContigs = referencePanelContigs,
			branch = branch,
			haplotype_database = haplotype_database
	}

	call ValidateScoring.ValidateScoring {
		input:
			validationArrays = validateImputation.imputed_multisample_vcf,
			validationArraysMain = validateImputation.imputed_multisample_vcf_main,
			validationWgs = validationWGS,
			population_basename = population_basename,
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			weights = weights,
			sample_name_map = sample_name_map,
			wgs_vcf_to_plink_mem = wgs_vcf_to_plink_mem,
			branch = branch
	}

	output {
		File correlations = validateImputation.correlations
		File accuracy = validateImputation.accuracy
		File correlations_plot = validateImputation.correlations_plot
		File aggregated_imputation_metrics = validateImputation.aggregated_imputation_metrics

		File? score_comparison_branch = ValidateScoring.score_comparison_branch
		File? score_comparison_main_vs_branch = ValidateScoring.score_comparison_main_vs_branch
		File? pc_plot = ValidateScoring.pc_plot
		File? raw_score_comparison_branch = ValidateScoring.raw_score_comparison_branch
		Int n_original_sites = ValidateScoring.n_original_sites
		Int n_subset_sites = ValidateScoring.n_subset_sites
		Int? n_subset_sites_wgs = ValidateScoring.n_subset_sites_wgs
	}
}

task GetSingleSampleName {
	input {
		File vcf
	}

	Int disk_size = 100 + ceil(size(vcf, "GB"))

	command <<<
		bcftools query -l ~{vcf} -o samples.list
	>>>

	runtime {
		docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
		memory: "3 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		String sampleName = read_lines("samples.list")[0]
	}
}

task SplitMultiSampleVcf {
	input {
		File multiSampleVcf
		String sample
	}

	parameter_meta {
		multiSampleVcf : {
			localization_optional : true
		}
	}

	command <<<
		gatk SelectVariants -V ~{multiSampleVcf} -sn ~{sample} -O ~{sample}.vcf.gz
	>>>

	runtime {
		docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File single_sample_vcf = "~{sample}.vcf.gz"
		File single_sample_vcf_index = "~{sample}.vcf.gz.tbi"
	}
}