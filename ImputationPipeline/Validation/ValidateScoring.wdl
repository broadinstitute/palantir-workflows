version 1.0

import "../ScoringPart.wdl" as Scoring

workflow ValidateScoring {
	input {
		Array[File] imputedArrays
		Array[File] wgssForScoring

		String population_basename

		File population_loadings
		File population_meansd
		File population_pcs
		File pruning_sites_for_pca # and the sites used for PCA
		File population_vcf

		File weights
		File sample_name_map
	}

	scatter (i in range(length(imputedArrays))) {
		File imputedArray = imputedArrays[i]
		File wgsForScoring = wgssForScoring[i]

		call ExtractIDs as extractImputedIDs {
			input:
				vcf = imputedArray,
				output_basename = "imputed"
		}

		call ExtractIDs as extractWGSIDs {
			input:
				vcf = wgsForScoring,
				output_basename = "wgs"
		}

		call SubsetWeights as SubsetWeightsImputed {
			input:
				sites = extractImputedIDs.ids,
				weights = weights
		}

		call SubsetWeights as SubsetWeightsWGS {
			input:
				sites = extractWGSIDs.ids,
				weights = SubsetWeightsImputed.subset_weights
		}

		call Scoring.ScoringImputedDataset as ScoreImputedSubset {
			input:
				weights = SubsetWeightsWGS.subset_weights,
				imputed_array_vcf = imputedArray,
				population_basename = population_basename,
				basename = "imputed",
				population_loadings = population_loadings,
				population_meansd = population_meansd,
				population_pcs = population_pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				population_vcf = population_vcf,
				redoPCA = true
		}

		call Scoring.ScoringImputedDataset as ScoreImputed {
			input:
				weights = weights,
				imputed_array_vcf = imputedArray,
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
				weights = SubsetWeightsWGS.subset_weights,
				imputed_array_vcf = wgsForScoring,
				population_basename = population_basename,
				basename = "imputed",
				population_loadings = population_loadings,
				population_meansd = population_meansd,
				population_pcs = population_pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				population_vcf = population_vcf,
				redoPCA = true
		}
	}

	call CompareScores as CompareScoresSubset {
		input:
			arrayScores = ScoreImputedSubset.adjusted_array_scores,
			wgsScores = ScoreWGS.adjusted_array_scores,
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
		Int n_original_sites = SubsetWeightsImputed.n_original_sites[0]
		Array[Int] n_subset_sites = SubsetWeightsImputed.n_subset_sites
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
		Array[File] arrayScores
		Array[File] wgsScores
		File sample_name_map
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(ggplot2)

		array_scores <- list("~{sep='","' arrayScores}") %>% map(read_tsv) %>% reduce(bind_rows) %>% transmute(IID, adjusted_score_array=adjusted_score)
		wgs_score <- list("~{sep='","' wgsScores}") %>% map(read_tsv) %>% reduce(bind_rows) %>% transmute(IID, adjusted_score_wgs=adjusted_score)

		sample_names <- read_delim("~{sample_name_map}", delim=":")

		combined_scores <- inner_join(inner_join(array_scores, sample_names, by=c("IID"="V1")), wgs_score, by=c("V2"="IID"))

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