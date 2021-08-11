version 1.0
import "PRSTasks.wdl" as PRSTasks
import "ScoringPart.wdl" as Score

workflow PRSWrapper {
  input {
    Array[String] condition_names
    Array[Boolean] use_condition
    Array[Float] percentile_thresholds
    Array[File] weights_files

    File vcf
    String sample_id

    File population_loadings
    File population_meansd
    File population_pcs
    File pruning_sites_for_pca # and the sites used for PCA
    File population_vcf
  }

  if (length(condition_names) != length(use_condition) || length(condition_names) != length(weights_files) || length(condition_names) != length(percentile_thresholds)) {
    call PRSTasks.ErrorWithMessage {
      input:
        message = "conditions_names, use_condition, use_ancestry_correction, and weights_files must all be arrays of the same length"
    }
  }

  scatter(i in range(length(condition_names))) {
    if (use_condition[i]) {
      call Score.ScoringImputedDataset {
        input:
          weights = weights_files[i],
          imputed_array_vcf = vcf,
          population_loadings = population_loadings,
          population_meansd = population_meansd,
          population_pcs = population_pcs,
          pruning_sites_for_pca = pruning_sites_for_pca,
          population_vcf = population_vcf
      }

      call SelectValuesOfInterest {
        input:
          score_result = ScoringImputedDataset.adjusted_array_scores,
          sample_id = sample_id,
          condition_name = condition_names[i],
          threshold = percentile_thresholds[i]
      }
    }
  }
}

task SelectValuesOfInterest {
  input {
    File score_result
    String sample_id
    String condition_name
    Float threshold
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    score <- read_tsv("~{score_result}")
    if (nrow(score) != 1) {
      quit(1)
    }
    if ((score %>% pull(`#IID`))[[1]] != ~{sample_id}) {
      quit(1)
    }

    raw_score <- (score %>% pull(SCORE1_SUM))[[1]]
    adjusted_score <- (score %>% pull(adjusted_score))[[1]]
    percentile <- (score %>% pull(percentile))[[1]]

    result <- tibble(sample_id = "~{sample_id}", ~{condition_name}_raw = raw_score, ~{condition_name}_adjusted = adjusted_score, ~{condition_name}_high = (percentile > threshold))
    write_csv(result, "result.tsv")

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "16 GB"
  }

  output {
    File raw_score_comparison_branch = "raw_score_comparison_~{branch}.png"
  }
}
