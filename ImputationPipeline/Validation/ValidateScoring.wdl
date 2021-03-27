version 1.0

import "../ScoringPart.wdl" as Scoring

workflow ValidateScoring {
	input {
		File validationArrays
		File validationWgs

		String population_basename

		File population_loadings
		File population_meansd
		File population_pcs
		File pruning_sites_for_pca # and the sites used for PCA
		File population_vcf

		File weights
		File sample_name_map
	}

	call ExtractIDs as extractImputedIDs {
		input:
			vcf = validationArrays,
			output_basename = "imputed"
	}


	call SubsetWeights {
		input:
			sites = extractImputedIDs.ids,
			weights = weights
	}


	call Scoring.ScoringImputedDataset as ScoreImputed {
		input:
			weights = weights,
			imputed_array_vcf = validationArrays,
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true
	}

	call Scoring.ScoringImputedDataset as ScoreWGSSubset {
		input:
			weights = SubsetWeights.subset_weights,
			imputed_array_vcf = validationWgs,
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true
	}

	call Scoring.ScoringImputedDataset as ScoreWGS {
		input:
			weights = weights,
			imputed_array_vcf = validationWgs,
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true
	}

	call CompareScores as CompareScoresSubset {
		input:
			arrayScores = ScoreImputed.adjusted_array_scores,
			wgsScores = ScoreWGSSubset.adjusted_array_scores,
			sample_name_map = sample_name_map
		}

	call CompareScores {
		input:
			arrayScores = ScoreImputed.adjusted_array_scores,
			wgsScores = ScoreWGS.adjusted_array_scores,
			sample_name_map = sample_name_map
		}


	output {
		File score_comparison_subset = CompareScoresSubset.score_comparison
		File score_comparison = CompareScores.score_comparison
		Int n_original_sites = SubsetWeights.n_original_sites
		Int n_subset_sites = SubsetWeights.n_subset_sites
	}
}


task SubsetWeights {
	input {
		File sites
		File weights
	}

	command <<<

	Rscript -<< "EOF"
	library(readr)
	library(dplyr)

	weights <- read_tsv("~{weights}")
	sites <- read_tsv("~{sites}", col_names = "variant")

	weights_subset <- inner_join(weights, sites)
	write_tsv(weights_subset, "weights_subset.txt")
	write_tsv(weights %>% count(), "n_original_weights.txt", col_names = FALSE)
	write_tsv(weights_subset %>% count(), "n_subset_weights.txt", col_names = FALSE)

	EOF

	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File subset_weights = "weights_subset.txt"
		Int n_original_sites = read_int("n_original_weights.txt")
		Int n_subset_sites = read_int("n_subset_weights.txt")
	}
}

task ExtractIDs {
	 input {
		 File vcf
		 String output_basename
		 Int disk_size = 2*ceil(size(vcf, "GB")) + 100
	 }

	 command <<<
		bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.original_array.ids
	 >>>
	 output {
		 File ids = "~{output_basename}.original_array.ids"
	 }
	 runtime {
		 docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
		 disks: "local-disk " + disk_size + " HDD"
		 memory: "4 GB"
	 }
}

task CompareScores {
	input {
		File arrayScores
		File wgsScores
		File sample_name_map
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(ggplot2)
		library(purrr)

		array_scores <- read_tsv("~{arrayScores}") %>% transmute(IID, adjusted_score_array=adjusted_score)
		wgs_score <- read_tsv("~{wgsScores}") %>% transmute(IID, adjusted_score_wgs=adjusted_score)

		sample_names <- read_delim("~{sample_name_map}", delim=":", col_names=FALSE)

		combined_scores <- inner_join(inner_join(array_scores, sample_names, by=c("IID"="X1")), wgs_score, by=c("X2"="IID"))

		ggplot(combined_scores, aes(x=adjusted_score_array, y=adjusted_score_wgs)) +
		geom_point() +
		geom_abline(intercept=0, slope=1) +
		xlab("Array Score") +
		ylab("WGS Score")

		ggsave(filename="score_comparison.png")

		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File score_comparison = "score_comparison.png"
	}
}