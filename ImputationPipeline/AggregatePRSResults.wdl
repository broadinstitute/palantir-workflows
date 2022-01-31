version 1.0

workflow AggregatePRSResults {
  input {
    Array[File] results
    Array[File] target_pc_projections
    File population_pc_projections
    String population_name = "Reference Population"
    File expected_control_results
  }

  call AggregateResults {
    input:
      results = results
  }

  call PlotPCA {
    input:
      batch_id = AggregateResults.batch_id,
      population_name = population_name,
      target_pc_projections = target_pc_projections,
      population_pc_projections = population_pc_projections
  }

  call BuildHTMLReport {
    input:
      batch_id = AggregateResults.batch_id,
      batch_all_results = AggregateResults.batch_all_results,
      batch_control_results = AggregateResults.batch_control_results,
      expected_control_results = expected_control_results,
      batch_summarised_results = AggregateResults.batch_summarised_results,
      score_distribution = AggregateResults.batch_score_distribution,
      pc_plot = PlotPCA.pc_plot
  }

  output {
    File batch_all_results = AggregateResults.batch_all_results
    File batch_control_results = AggregateResults.batch_control_results
    File batch_summarised_results = AggregateResults.batch_summarised_results
    File score_distribution = AggregateResults.batch_score_distribution
    File pc_polot = PlotPCA.pc_plot
    File report = BuildHTMLReport.report
  }
}

task AggregateResults {
  input {
    Array[File] results
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
    
    num_batches <- lapply(results, function(row) {length(unique(row))})$batch_id

    if (num_batches != 1) {
      stop(paste0("There are ", num_batches, " batch IDs in the input tables, however, only 1 is expected."))
    }

    batch_id <- results$batch_id[1]

    num_control_samples <- sum(results$is_control=="true")

    if (num_batches != 1) {
      stop(paste("There are ", num_control_samples, " control samples in the input tables, however, only 1 is expected."))
    }

    write_tsv(results, paste0(batch_id, "_all_results.tsv"))

    write_tsv(results %>% filter(is_control == "true"), paste0(batch_id, "_control_results.tsv"))

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

    write_tsv(results_summarised, paste0(batch_id, "_summarised_results.tsv"))

    ggplot(results_pivoted, aes(x=adjusted)) +
      geom_density(aes(color=condition), fill=NA, position = "identity") +
      xlim(-10,10) + theme_bw() + xlab("z-score") + geom_function(fun=dnorm) +
      ylab("density")
    ggsave(filename = paste0(batch_id, "_score_distribution.png"), dpi=300, width = 6, height = 6)

    writeLines(batch_id, "batch_id.txt")

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    String batch_id = read_string("batch_id.txt")
    File batch_all_results = glob("*_all_results.tsv")[0]
    File batch_control_results = glob("*_control_results.tsv")[0]
    File batch_summarised_results = glob("*_summarised_results.tsv")[0]
    File batch_score_distribution = glob("*_score_distribution.png")[0]
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

    ggplot(population_pcs, aes(x=PC1, y=PC2, color="~{population_name}")) +
      geom_point(size=0.1, alpha=0.1) +
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

task BuildHTMLReport {
  input {
    File batch_all_results
    File batch_control_results
    File expected_control_results
    File batch_summarised_results
    File score_distribution
    File pc_plot
    String batch_id
  }

  command <<<
    set -xeo pipefail

    cat << EOF > ~{batch_id}_report.Rmd
    ---
    title: "Batch ~{batch_id} PRS Summary"
    output:
    html_document:
      df_print: paged
    date: "$(date)"
    ---

    \`\`\`{r setup, include=FALSE}
    library(readr)
    library(ggplot2)
    library(knitr)
    library(dplyr)
    library(stringr)
    library(purrr)
    batch_all_results <- read_tsv("~{batch_all_results}")
    batch_control_results <- read_tsv("~{batch_control_results}")
    expected_control_results <- read_tsv("~{expected_control_results}")
    batch_summary <- read_tsv("~{batch_summarised_results}")
    batch_summary <- batch_summary %>% rename_with(.cols = -condition, ~ str_to_title(gsub("_"," ", .x)))
    \`\`\`


    ## Control Sample
    ```{r control, echo = FALSE, results = "asis"}
    kable(list(batch_control_results, expected_control_results) %>% reduce(bind_rows) %>% select(-batch_id) %>% select(-is_control) , digits = 2)
    ```

    ## Batch Summary
    \`\`\`{r summary table, echo = FALSE, results = "asis" }
    kable(batch_summary, digits = 2)
    \`\`\`



    ## Batch Score distribution
    ![](~{score_distribution})

    ## PCA
    ![](~{pc_plot})

    ## Individual Sample Results
    \`\`\`{r sample results , echo = FALSE, results = "asis"}
    kable(batch_all_results %>% mutate(across(ends_with("risk"), ~ kableExtra::cell_spec(.x, color=ifelse(is.na(.x), "blue", ifelse(.x=="NOT_RESULTED", "red", ifelse(.x == "HIGH", "orange", "green")))))), digits = 2)
    \`\`\`
    EOF

    Rscript -e "library(rmarkdown); rmarkdown::render('~{batch_id}_report.Rmd', 'html_document')"
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/tidyverse_kableextra_docker@sha256:e00dd57a1c446e9c07429929893a8399c4a299f0b4fd182245bf45e77c6ab8bd"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File report = "~{batch_id}_report.html"
  }
}
