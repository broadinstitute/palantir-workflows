version 1.0
import "ScoringPart.wdl" as Score
import "CKDRiskAdjustment.wdl" as CKDRiskAdjustmentWF
import "Structs.wdl"

workflow PRSWrapper {
  input {
    Array[PRSWrapperConditionResource] condition_resources
    File? ckd_risk_alleles
    Boolean ckd_adjust_monogenic_risk = true
    Float z_score_reportable_range

    File vcf
    String sample_set_id
    String lab_batch_id
    Boolean is_control_sample_in
    Boolean redoPCA = false

    File? population_loadings
    File? population_meansd
    File? population_pcs
    File? pruning_sites_for_pca # and the sites used for PCA
    Boolean adjust_scores = true
    Int mem_extract
    Int mem_vcf_to_plink

    Boolean use_ref_alt_for_ids
  }

  scatter(condition_resource in condition_resources) {
    if (condition_resource.score_condition) {
      call Score.ScoringImputedDataset {
        input:
          named_weight_set = condition_resource.named_weight_set,
          imputed_array_vcf = vcf,
          adjustScores = adjust_scores,
          population_loadings = population_loadings,
          population_meansd = population_meansd,
          pruning_sites_for_pca = pruning_sites_for_pca,
          population_pcs = population_pcs,
          basename = sample_set_id,
          fitted_model_params_and_sites = condition_resource.ancestry_model_params_and_sites,
          redoPCA = redoPCA,
          mem_extract = mem_extract,
          vcf_to_plink_mem = mem_vcf_to_plink,
          use_ref_alt_for_ids = use_ref_alt_for_ids
      }
    }
  }


  output {
    Array[File?] raw_scores = select_all(ScoringImputedDataset.raw_scores)
  }
}


task CheckZScoreAgainstReportableRange {
  input {
    File score_result
    Float z_score_reportable_range
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    score <- read_tsv("~{score_result}")
    if (nrow(score) != 1) {
    quit(status=1)
    }

    adjusted_score <- (score %>% pull(adjusted_score))[[1]]
    write(abs(adjusted_score) > ~{z_score_reportable_range}, "out_of_reportable_range.bool")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    Boolean out_of_reportable_range = read_boolean("out_of_reportable_range.bool")
  }
}

task SelectValuesOfInterest {
  input {
    File score_result
    String sample_id
    String condition_name
    Float threshold
    Boolean out_of_reportable_range
    Float z_score_reportable_range
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
    risk <- ifelse(percentile > ~{threshold}, "HIGH", "NOT_HIGH")

    raw_score_output <- ~{if out_of_reportable_range then '"NOT_RESULTED"' else 'raw_score'}
    adjusted_score_output <- ~{if out_of_reportable_range then '"NOT_RESULTED"' else 'adjusted_score'}
    percentile_output <- ~{if out_of_reportable_range then '"NOT_RESULTED"' else 'percentile'}
    risk_output <- ~{if out_of_reportable_range then '"NOT_RESULTED"' else 'risk'}
    reason_not_resulted <- ~{if out_of_reportable_range then
                                'ifelse(adjusted_score > 0, "Z-SCORE ABOVE + ' + z_score_reportable_range +'", "Z-SCORE BELOW - ' + z_score_reportable_range +'")' else
                                "NA"
                           }

    result <- tibble(sample_id = "~{sample_id}", ~{condition_name}_raw = raw_score_output,
                                                 ~{condition_name}_adjusted = adjusted_score_output,
                                                 ~{condition_name}_percentile = percentile_output,
                                                 ~{condition_name}_risk = risk_output,
                                                 ~{condition_name}_reason_not_resulted = reason_not_resulted)
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
    echo "sample_id, ~{condition_name}_raw, ~{condition_name}_adjusted, ~{condition_name}_percentile, ~{condition_name}_risk, ~{condition_name}_reason_not_resulted" > results.csv
    echo "~{sample_id}, NA, NA, NA, NA, NA" >> results.csv
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
    String lab_batch
    Boolean is_control_sample
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(tibble)

    results <- c("~{sep = '","' results_in}") %>% map(read_csv) %>% reduce(inner_join)
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

task CombineMissingSitesShiftedScores {
  input {
    Array[File] missing_sites_shifted
    String lab_batch
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)

    missing_sites_shifted <- c("~{sep = '","' missing_sites_shifted}") %>% map(read_tsv) %>% reduce(bind_rows) %>% rename(sample_id = IID)
    write_tsv(missing_sites_shifted, "~{lab_batch}.missing_sites_shifts.tsv")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File missing_sites_shifts = "~{lab_batch}.missing_sites_shifts.tsv"
  }
}
