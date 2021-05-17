version 1.0

import "../ScoringPart.wdl" as Scoring
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/ImputationPipeline/ScoringPart.wdl" as ScoringMain

workflow ValidateScoring {
	input {
		File? validationArrays #array to score with this branch
		File validationArraysMain #array to score with main branch (will be used for this branch also if validationArrays not provided)
		File validationWgs #wgs to score for comparison

		String population_basename #sets name of population related output files

		File population_loadings #population pca
		File population_meansd #population pca
		File population_pcs #population pca
		File pruning_sites_for_pca # and the sites used for PCA
		File population_vcf

		File weights
		File sample_name_map #file mapping sample names in arrays to sample names in wgs.  ":" used as separator
		String branch #name of branch being tested.  used for display in plots, does not effect computation

		Int wgs_vcf_to_plink_mem = 8
		Boolean adjustScores = true
	}

	#Extract the sites which are in the imputed vcf.  We will score the wgs over only the sites included in the imputed vcf
	call ExtractIDs as extractImputedIDs {
		input:
			vcf = select_first([validationArrays, validationArraysMain]),
			output_basename = "imputed"
	}


	#Subset the weights to only the sites in the imputed vcf
	call SubsetWeights {
		input:
			sites = extractImputedIDs.ids,
			weights = weights
	}

	#remove sites with any no-call genotypes from the wgs scoring.  otherwise no-calls in the wgs would skew the wgs scores in a way that would not be accounted for by the score adjustment
	call QCSites {
		input:
			input_vcf = validationWgs,
			output_vcf_basename = "wgsValidation"
	}

	if (!adjustScores) {
		#if not adjusting scores, must subset weights also to only sites with no-calls in wgs
		call ExtractIDs as extractWGSIDs {
			input:
				vcf = QCSites.output_vcf,
				output_basename = "wgs"
		}

		call SubsetWeights as SubsetWeightsWGS {
			input:
				sites = extractWGSIDs.ids,
				weights = SubsetWeights.subset_weights
		}
	}

	#run scoring on this branch, using imputed data from this branch, or shared imputed data is we are studying only changes in scoring
	call Scoring.ScoringImputedDataset as ScoreImputed {
		input:
			weights = select_first([SubsetWeightsWGS.subset_weights, weights]),
			imputed_array_vcf = select_first([validationArrays, validationArraysMain]),
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true,
			adjustScores = adjustScores
	}

	#run scoring on main branch
	call ScoringMain.ScoringImputedDataset as ScoreImputedMain {
		input:
			weights = weights,
			imputed_array_vcf = validationArraysMain,
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true
	}

	#score wgs over only sites in the imputed array which are called in every wgs sample
	call Scoring.ScoringImputedDataset as ScoreWGS {
		input:
			weights = select_first([SubsetWeightsWGS.subset_weights, SubsetWeights.subset_weights]),
			imputed_array_vcf = QCSites.output_vcf,
			population_basename = population_basename,
			basename = "imputed",
			population_loadings = population_loadings,
			population_meansd = population_meansd,
			population_pcs = population_pcs,
			pruning_sites_for_pca = pruning_sites_for_pca,
			population_vcf = population_vcf,
			redoPCA = true,
			vcf_to_plink_mem = wgs_vcf_to_plink_mem,
			adjustScores = adjustScores
	}

	if (adjustScores) {
		#compare this branch scores to wgs scores and to main branch scores
		call CompareScores {
			input:
				arrayScores = select_first([ScoreImputed.adjusted_array_scores]),
				wgsScores = select_first([ScoreWGS.adjusted_array_scores]),
				arrayScoresMain = ScoreImputedMain.adjusted_array_scores,
				branch = branch,
				sample_name_map = sample_name_map
			}
	}

	if (!adjustScores) {
		#compare raw scores to wgs.  will add in comparison to main branch later
		call CompareRawScores {
			input:
				arrayScores = ScoreImputed.raw_scores,
				wgsScores = ScoreWGS.raw_scores,
				branch = branch,
				sample_name_map = sample_name_map
		}
	}


	output {
		File? score_comparison_branch = CompareScores.score_comparison_branch
		File? score_comparison_main_vs_branch = CompareScores.score_comparison_main_vs_branch

		File? pc_plot = ScoreImputed.pc_plot
		File? raw_score_comparison_branch = CompareRawScores.raw_score_comparison_branch
		Int n_original_sites = SubsetWeights.n_original_sites
		Int n_subset_sites = SubsetWeights.n_subset_sites
		Int? n_subset_sites_wgs = SubsetWeightsWGS.n_subset_sites
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
		File arrayScoresMain
		File wgsScores
		File sample_name_map
		String branch
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(ggplot2)
		library(purrr)

		array_scores <- read_tsv("~{arrayScores}") %>% transmute(IID, adjusted_score_array=adjusted_score)
		wgs_score <- read_tsv("~{wgsScores}") %>% transmute(IID, adjusted_score_wgs=adjusted_score)

		array_scores_main <- read_tsv("~{arrayScoresMain}") %>% transmute(IID, adjusted_score_array=adjusted_score)

		sample_names <- read_delim("~{sample_name_map}", delim=":", col_names=FALSE)



		combined_scores <- inner_join(inner_join(array_scores, sample_names, by=c("IID"="X1")), wgs_score, by=c("X2"="IID"))

		score_main_branch <- inner_join(array_scores %>% transmute(IID, adjusted_score_branch=adjusted_score_array), array_scores_main %>% transmute(IID, adjusted_score_main=adjusted_score_array))

		ggplot(combined_scores, aes(x=adjusted_score_array, y=adjusted_score_wgs)) +
		geom_point() +
		geom_abline(intercept=0, slope=1) +
		xlab("Array Score ~{branch}") +
		ylab("WGS Score ~{branch}")

		ggsave(filename="score_comparison_~{branch}.png")

		ggplot(score_main_branch, aes(x=adjusted_score_branch, y=adjusted_score_main)) +
		geom_point() +
		geom_abline(intercept=0, slope=1) +
		xlab("Array Score ~{branch}") +
		ylab("Array Score Main")

		ggsave(filename="score_comparison_main_vs_branch.png")

		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File score_comparison_branch = "score_comparison_~{branch}.png"
		File score_comparison_main_vs_branch = "score_comparison_main_vs_branch.png"
	}
}

task CompareRawScores {
	input {
		File arrayScores
		File wgsScores
		File sample_name_map
		String branch
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(ggplot2)
		library(purrr)

		array_scores <- read_tsv("~{arrayScores}") %>% transmute(`#IID`, raw_score_array=SCORE1_SUM)
		wgs_score <- read_tsv("~{wgsScores}") %>% transmute(`#IID`, raw_score_wgs=SCORE1_SUM)

		sample_names <- read_delim("~{sample_name_map}", delim=":", col_names=FALSE)



		combined_scores <- inner_join(inner_join(array_scores, sample_names, by=c("IID"="X1")), wgs_score, by=c("X2"="IID"))

		ggplot(combined_scores, aes(x=raw_score_array, y=raw_score_wgs)) +
		geom_point() +
		geom_abline(intercept=0, slope=1) +
		xlab("Raw Array Score ~{branch}") +
		ylab("Raw WGS Score ~{branch}")

		ggsave(filename="raw_score_comparison_~{branch}.png")

		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "16 GB"
	}

	output {
		File raw_score_comparison_branch = "raw_score_comparison_~{branch}.png"
	}
}

task QCSites {
  input {
    File input_vcf
    String output_vcf_basename
   }
    Int disk_size = ceil(2*(size(input_vcf, "GB")))

  command <<<
    # site with any missing genotypes are removed;
    vcftools --gzvcf ~{input_vcf}  --max-missing-count 0 --recode -c | bgzip -c > ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "skwalker/imputation:with_vcftools" # TODO: use a public one (not suse bcftools biocontainers also has vcftools )
    memory: "16 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}