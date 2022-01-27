version 1.0
import "ScoringPart.wdl" as Score
import "CKDRiskAdjustment.wdl" as CKDRiskAdjustment

workflow PRSWrapper {
  input {
    Array[String] condition_names
    Array[Boolean] score_condition
    Array[Float] percentile_thresholds
    Array[File] weights_files
    File ckd_risk_alleles

    File vcf
    String sample_id
    String batch_id
    Boolean is_control
    Boolean redoPCA = false

    File population_loadings
    File population_meansd
    File population_pcs
    File pruning_sites_for_pca # and the sites used for PCA
    File population_vcf
    String population_basename
  }

  if (length(condition_names) != length(score_condition) || length(condition_names) != length(weights_files) || length(condition_names) != length(percentile_thresholds)) {
    call ErrorWithMessage {
      input:
        message = "conditions_names, use_condition, use_ancestry_correction, and weights_files must all be arrays of the same length"
    }
  }

  scatter(i in range(length(condition_names))) {
    if (score_condition[i]) {
      call Score.ScoringImputedDataset {
        input:
          weights = weights_files[i],
          imputed_array_vcf = vcf,
          population_loadings = population_loadings,
          population_meansd = population_meansd,
          population_pcs = population_pcs,
          pruning_sites_for_pca = pruning_sites_for_pca,
          population_vcf = population_vcf,
          basename = sample_id,
          population_basename = population_basename,
          redoPCA = redoPCA
      }

      if (condition_names[i] == "ckd") {
        call CKDRiskAdjustment.CKDRiskAdjustment {
          input:
            adjustedScores = select_first([ScoringImputedDataset.adjusted_array_scores]),
            vcf = vcf,
            risk_alleles = ckd_risk_alleles
        }
      }


      call SelectValuesOfInterest {
        input:
          score_result = select_first([CKDRiskAdjustment.adjusted_scores_with_apol1, ScoringImputedDataset.adjusted_array_scores]),
          sample_id = sample_id,
          condition_name = condition_names[i],
          threshold = percentile_thresholds[i]
      }
    }

    if (!score_condition[i])
    {
      call CreateUnscoredResult {
        input:
          sample_id = sample_id,
          condition_name = condition_names[i]
      }
    }

    File result_for_condition = select_first([SelectValuesOfInterest.results, CreateUnscoredResult.results])
  }

  call JoinResults{
    input:
      results_in = result_for_condition,
      batch_id = batch_id,
      is_control = is_control
  }

  output {
    File results = JoinResults.results
    File pcs = select_first(ScoringImputedDataset.pc_projection)
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
      quit(status=1)
    }
    if ((score %>% pull(`IID`))[[1]] != "~{sample_id}") {
      quit(status=1)
    }

    raw_score <- (score %>% pull(SCORE1_SUM))[[1]]
    adjusted_score <- (score %>% pull(adjusted_score))[[1]]
    percentile <- (score %>% pull(percentile))[[1]]

    result <- tibble(sample_id = "~{sample_id}", ~{condition_name}_raw = raw_score, ~{condition_name}_adjusted = adjusted_score, ~{condition_name}_percentile = percentile, ~{condition_name}_risk = ifelse(percentile > ~{threshold}, "HIGH", "NOT_HIGH"))
    write_csv(result, "results.csv")

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File results = "results.csv"
  }
}

task CreateUnscoredResult {
  input {
    String condition_name
    String sample_id
  }

  command <<<
    echo "sample_id, ~{condition_name}_raw, ~{condition_name}_adjusted, ~{condition_name}_percentile", ~{condition_name}_risk > results.csv
    echo "~{sample_id}, NA, NA, NA, NA" >> results.csv
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File results = "results.csv"
  }
}

task JoinResults {
  input {
    Array[File] results_in
    String batch_id
    Boolean is_control
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)

    results <- c("~{sep = '","' results_in}") %>% map(read_csv) %>% reduce(inner_join) %>% add_column(batch_id=~{batch_id}, .after="sample_id") %>% add_column(batch_id=~{is_control}, .after="batch_id")
    write_csv(results, "results.csv")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File results = "results.csv"
  }
}

#Print given message to stderr and return an error
task ErrorWithMessage{
  input {
    String message
  }
  command <<<
    >&2 echo "Error: ~{message}"
    exit 1
  >>>

  runtime {
    docker: "ubuntu"
  }
}
