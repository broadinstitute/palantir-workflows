version 1.0

workflow AggregatePRSResults {
  input {
    String batch_id
    Array[File] results
    Array[File] target_pc_projections
    File population_pc_projections
    String population_name = "Reference Population"
  }

  call AggregateResults {
    input:
      results = results,
      batch_id = batch_id
  }

  call PlotPCA {
    input:
      batch_id = batch_id,
      population_name = population_name,
      target_pc_projections = target_pc_projections,
      population_pc_projections = population_pc_projections
  }

  output {
    File batch_results = AggregateResults.batch_results
    File batch_summarised_results = AggregateResults.batch_summarised_results
    File score_distribution = AggregateResults.score_distribution
    File pc_polot = PlotPCA.pc_plot
  }
}

task AggregateResults {
  input {
    Array[File] results
    String batch_id
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(magrittr)
    library(ggplot2)

    results <- c("~{sep='","' results}") %>% map(read_csv, col_types=cols(.default = 'c')) %>% reduce(bind_rows)
    write_tsv(results, "~{batch_id}_results.tsv")

    results_pivoted <- results %>% pivot_longer(!sample_id, names_to=c("condition",".value"), names_pattern="([^_]+)_(.+)")
    results_pivoted <- results_pivoted %T>% {options(warn=-1)} %>% mutate(adjusted = as.numeric(adjusted),
                                                                          raw = as.numeric(raw),
                                                                          percentile = as.numeric(percentile)) %T>% {options(warn=0)}

    results_summarised <- results_pivoted %>% group_by(condition) %>%
                                              summarise(across(c(adjusted,percentile), ~mean(.x, na.rm=TRUE), .names = "mean_{.col}"),
                                                        num_samples=n(),
                                                        num_scored = sum(!is.na(risk)),
                                                        num_high = sum(risk=="HIGH", na.rm=TRUE),
                                                        num_not_high = sum(risk=="NOT_HIGH", na.rm=TRUE),
                                                        num_not_resulted = sum(risk=="NOT_RESULTED", na.rm = TRUE))

    write_tsv(results_summarised, "~{batch_id}_summarised_results.tsv")

    ggplot(results_pivoted, aes(x=adjusted)) +
      geom_density(aes(color=condition), fill=NA, position = "identity") +
      xlim(-10,10) + theme_bw() + xlab("z-score") + geom_function(fun=dnorm) +
      ylab("density")
    ggsave(filename = "~{batch_id}_score_distribution.png", dpi=300, width = 6, height = 6)

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File batch_results = "~{batch_id}_results.tsv"
    File batch_summarised_results = "~{batch_id}_summarised_results.tsv"
    File score_distribution = "~{batch_id}_score_distribution.png"
  }
}

task PlotPCA {
  input {
    Array[File] target_pc_projections
    File population_pc_projections
    String batch_id
    String population_name
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(ggplot2)

    target_pcs <- c("~{sep='","' target_pc_projections}") %>% map(read_tsv) %>% reduce(bind_rows)
    population_pcs <- read_tsv("~{population_pc_projections}")

    ggplot(population_pcs, aes(x=PC1, y=PC2, color="~{population_name}"), size=0.1, alpha=0.1) +
      geom_point() +
      geom_point(data=target_pcs, aes(color="~{batch_id}")) +
      theme_bw()

    ggsave(filename = "~{batch_id}_PCA_plot.png", dpi=300, width = 6, height = 6)

    EOF

  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File pc_plot = "~{batch_id}_PCA_plot.png"
  }
}
