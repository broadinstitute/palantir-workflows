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
    }
  }
}

task SelectValuesOfInterest {
  input {
    File score_result
    String sample_id
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    score <- read_tsv("~{score_result}")
    if (score %>% pull(`#IID`) == ~{sample_id}) {

    }
    EOF
  >>>
}
