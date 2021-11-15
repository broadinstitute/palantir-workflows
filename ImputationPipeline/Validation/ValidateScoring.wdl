version 1.0

import "../ScoringPart.wdl" as Scoring
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/ImputationPipeline/ScoringPart.wdl" as ScoringMain
import "../TrainAncestryAdjustmentModel.wdl" as TrainModel
import "../PCATasks.wdl" as PCATasks
import "../Structs.wdl" as Structs
import "SubsetWeightSet.wdl" as SubsetWeightSet

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

		Array[NamedWeightSet] named_weight_sets
		File sample_name_map #file mapping sample names in arrays to sample names in wgs.  ":" used as separator
		String branch #name of branch being tested.  used for display in plots, does not effect computation

		Int wgs_vcf_to_plink_mem = 8

		Float max_diff_pretrain_threshold = 0.0000000001 # maximum difference allowed between pretraining model separately and as part of scoring
	}

	#Extract the sites which are in the imputed vcf.  We will score the wgs over only the sites included in the imputed vcf
	call ExtractIDs as extractImputedIDs {
		input:
			vcf = select_first([validationArrays, validationArraysMain]),
			output_basename = "imputed"
	}

	#need to redo PCA in case input PCA not appropriate for this set of sites
	call PCATasks.ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
		input:
			vcf = select_first([population_vcf]),
			pruning_sites = select_first([pruning_sites_for_pca]),
			subset_to_sites = extractImputedIDs.ids,
			basename = "population"
	}

	call PCATasks.PerformPCA {
		input:
			bim = PopulationArrayVcfToPlinkDataset.bim,
			bed = PopulationArrayVcfToPlinkDataset.bed,
			fam = PopulationArrayVcfToPlinkDataset.fam,
			basename = "population"
	}

	scatter (named_weight_set in named_weight_sets) {
		String conditions = named_weight_set.condition_name
		#Subset the weights to only the sites in the imputed vcf
		call SubsetWeightSet.SubsetWeightSet {
			input:
				sites_to_subset_to = extractImputedIDs.ids,
				weight_set = named_weight_set.weight_set
		}

		#remove sites with any no-call genotypes from the wgs scoring.  otherwise no-calls in the wgs would skew the wgs scores in a way that would not be accounted for by the score adjustment
		call QCSites {
			input:
				input_vcf = validationWgs,
				output_vcf_basename = "wgsValidation"
		}

		#subset weights also to only sites without no-calls in wgs
		call ExtractIDs as extractWGSIDs {
			input:
				vcf = QCSites.output_vcf,
				output_basename = "wgs"
		}

		call SubsetWeightSet.SubsetWeightSet as SubsetWeightSetWGS {
			input:
				sites_to_subset_to = extractWGSIDs.ids,
				weight_set = SubsetWeightSet.subsetted_weight_set
		}

		#run scoring on this branch, using imputed data from this branch, or shared imputed data is we are studying only changes in scoring
		call Scoring.ScoringImputedDataset as ScoreImputed {
			input:
				named_weight_set = named_weight_set,
				imputed_array_vcf = select_first([validationArrays, validationArraysMain]),
				population_basename = population_basename,
				basename = "imputed",
				population_loadings = population_loadings,
				population_meansd = population_meansd,
				population_pcs = population_pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				population_vcf = population_vcf,
				redoPCA = true
		}

		call Scoring.ScoringImputedDataset as ScoreImputedForWGSComparison {
			input:
				named_weight_set = object{condition_name : conditions, weight_set : SubsetWeightSetWGS.subsetted_weight_set},
				imputed_array_vcf = select_first([validationArrays, validationArraysMain]),
				population_basename = population_basename,
				basename = "imputed",
				population_loadings = population_loadings,
				population_meansd = population_meansd,
				population_pcs = population_pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				population_vcf = population_vcf,
				redoPCA = true
		}

		#run scoring on main branch
		call ScoringMain.ScoringImputedDataset as ScoreImputedMain {
			input:
				weight_set = named_weight_set.weight_set,
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
				named_weight_set = object{condition_name : conditions, weight_set : SubsetWeightSetWGS.subsetted_weight_set},
				imputed_array_vcf = QCSites.output_vcf,
				population_basename = population_basename,
				basename = "imputed",
				population_loadings = population_loadings,
				population_meansd = population_meansd,
				population_pcs = population_pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				population_vcf = population_vcf,
				redoPCA = true,
				vcf_to_plink_mem = wgs_vcf_to_plink_mem
		}

		#train ancestry adjustment separately and score.  scores should be indentical to standard run

		call TrainModel.TrainAncestryAdjustmentModel {
			input:
				named_weight_set = named_weight_set,
				population_pcs = PerformPCA.pcs,
				population_vcf = population_vcf,
				population_basename = population_basename
		}

		call Scoring.ScoringImputedDataset as ScoreWithTrainingResults {
			input:
				named_weight_set = named_weight_set,
				imputed_array_vcf = validationArraysMain,
				basename = "imputed",
				population_loadings = PerformPCA.pc_loadings,
				population_meansd = PerformPCA.mean_sd,
				population_pcs = PerformPCA.pcs,
				pruning_sites_for_pca = pruning_sites_for_pca,
				fitted_model_params_and_sites = object{fitted_model_params : TrainAncestryAdjustmentModel.fitted_params, sites_used_in_scoring : TrainAncestryAdjustmentModel.sites_used_in_scoring}
		}
	}

	#compare this branch scores to wgs scores and to main branch scores
	call CompareScores as CompareScoresImputedWGS{
		input:
			x_axis_scores = select_all(ScoreImputedForWGSComparison.adjusted_array_scores),
			y_axis_scores = select_all(ScoreWGS.adjusted_array_scores),
			conditions = conditions,
			x_axis_label = "Array Score " + branch,
			y_axis_label = "WGS Score " + branch,
			sample_name_map = sample_name_map,
			output_filename = "score_comparison_" + branch
	}

	call CompareScores as CompareScoresMainVsBranch {
		input:
			x_axis_scores = select_all(ScoreImputed.adjusted_array_scores),
			y_axis_scores = select_all(ScoreImputedMain.adjusted_array_scores),
			conditions = conditions,
			x_axis_label = "Array Score " + branch,
			y_axis_label = "Array Score Main",
			output_filename = "score_comparison_main_vs_" + branch
	}

		call CompareScores as CompareScoresPreTraining {
			input:
				x_axis_scores = select_all(ScoreImputed.adjusted_array_scores),
				y_axis_scores = select_all(ScoreWithTrainingResults.adjusted_array_scores),
				conditions = conditions,
				x_axis_label = "Array Score " + branch,
				y_axis_label = "Array Score Pre-train model " + branch,
				output_filename = "score_comparison_pre_train_model_" + branch
		}

		if (CompareScoresPreTraining.max_score_diff > max_diff_pretrain_threshold) {
			call Scoring.ErrorWithMessage {
					input:
						message = "pre-trained score difference too large (" + CompareScoresPreTraining.max_score_diff + ")"
			}
		}

	output {
		File score_comparison_branch = CompareScoresImputedWGS.score_comparison
		File score_comparison_main_vs_branch = CompareScoresMainVsBranch.score_comparison
		File score_comparison_pretrian = CompareScoresPreTraining.score_comparison
		Float pretrain_max_diff = CompareScoresPreTraining.max_score_diff
		Float branch_main_max_diff = CompareScoresMainVsBranch.max_score_diff

		File pc_plot = select_first(ScoreImputed.pc_plot)
		Array[Int] n_original_weights = SubsetWeightSet.n_original_weights
		Array[Int] n_subset_weights = SubsetWeightSet.n_subset_weights
		Array[Int] n_subset_weights_wgs = SubsetWeightSetWGS.n_subset_weights
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
		Array[File] x_axis_scores
		Array[File] y_axis_scores
		Array[String] conditions
		File? sample_name_map
		String x_axis_label
		String y_axis_label
		String output_filename
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(ggplot2)
		library(purrr)

		conditions <- c("~{sep = '", "' conditions}")
		x_scores <- list("~{sep = '", "' x_axis_scores}") %>% map(read_tsv) %>% map2(conditions, ~ .x %>% mutate(condition = .y)) %>% reduce(bind_rows) %>% transmute(IID, condition, adjusted_score_x=adjusted_score)
		y_scores <- list("~{sep = '", "' y_axis_scores}") %>% map(read_tsv) %>% map2(conditions, ~ .x %>% mutate(condition = .y)) %>% reduce(bind_rows) %>% transmute(IID, condition, adjusted_score_y=adjusted_score)

		sample_names_map <- if(~{true="TRUE" false="FALSE" defined(sample_name_map)}) read_delim("~{sample_name_map}", delim=":", col_names=FALSE) else x_scores %>% transmute(X1=IID, X2=IID)

		combined_scores <- inner_join(inner_join(x_scores, sample_names_map, by=c("IID"="X1")), y_scores, by=c("X2"="IID", "condition"))

		ggplot(combined_scores, aes(x=adjusted_score_x, y=adjusted_score_y)) +
		geom_point() +
		geom_abline(intercept=0, slope=1) +
		facet_wrap(~condition) +
		xlab("~{x_axis_label}") +
		ylab("~{y_axis_label}")

		ggsave(filename="~{output_filename}.png")

		max_diff_conditions <- combined_scores %>% group_by(condition) %>% summarise(max_diff = max(abs(adjusted_score_x - adjusted_score_y)))
		write_tsv(max_diff_conditions, "max_diff_conditions.txt")
		max_diff <- max_diff_conditions %>% summarise(max_diff = max(max_diff))
		write_tsv(max_diff, "max_diff.txt", col_names = FALSE)
		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "4 GB"
	}

	output {
		File score_comparison = "~{output_filename}.png"
		Float max_score_diff = read_float("max_diff.txt")
		File max_diff_conditions = "max_diff_conditions.txt"
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
