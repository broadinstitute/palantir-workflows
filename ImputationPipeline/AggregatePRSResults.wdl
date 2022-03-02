version 1.0

workflow AggregatePRSResults {
  input {
    Array[File] results
    Array[File] target_pc_projections
    Array[File] missing_sites_shifts
    File population_pc_projections
    String population_name = "Reference Population"
    File expected_control_results
    String lab_batch
  }

  call AggregateResults {
    input:
      results = results,
      missing_sites_shifts = missing_sites_shifts,
      lab_batch = lab_batch
  }

  call PlotPCA {
    input:
      lab_batch = lab_batch,
      population_name = population_name,
      target_pc_projections = target_pc_projections,
      population_pc_projections = population_pc_projections
  }

  call BuildHTMLReport {
    input:
      lab_batch = lab_batch,
      batch_all_results = AggregateResults.batch_all_results,
      batch_control_results = AggregateResults.batch_control_results,
      batch_missing_sites_shifts = AggregateResults.batch_missing_sites_shifts,
      expected_control_results = expected_control_results,
      batch_summarised_results = AggregateResults.batch_summarised_results,
      score_distribution = AggregateResults.batch_score_distribution,
      target_pc_projections = target_pc_projections,
      population_pc_projections = population_pc_projections,
      population_name = population_name
  }

  output {
    File batch_all_results = AggregateResults.batch_all_results
    File batch_control_results = AggregateResults.batch_control_results
    File batch_summarised_results = AggregateResults.batch_summarised_results
    File batch_missing_sites_shifts =  AggregateResults.batch_missing_sites_shifts
    File score_distribution = AggregateResults.batch_score_distribution
    File pc_plot = PlotPCA.pc_plot
    File report = BuildHTMLReport.report
  }
}

task AggregateResults {
  input {
    Array[File] results
    Array[File] missing_sites_shifts
    String lab_batch
  }

  command <<<
    Rscript - <<- "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(magrittr)
    library(ggplot2)

    results <- c("~{sep='","' results}") %>% map(read_csv, col_types=cols(is_control_sample='l', .default='c')) %>% reduce(bind_rows)
    
    lab_batch <- results %>% pull(lab_batch) %>% unique()

    if (length(lab_batch) != 1) {
      stop(paste0("There are ", length(lab_batch), " lab batch IDs in the input tables, however, only 1 is expected."))
    }
    # If there is only one lab_batch, then lab_batch will have length 1 and can be treated as a single value from here on

    if (lab_batch != "~{lab_batch}") {
      stop(paste0("Expected lab batch ~{lab_batch} but found lab batch ", lab_batch))
    }

    num_control_samples <- results %>% filter(is_control_sample) %>% count()

    if (num_control_samples != 1) {
      stop(paste0("There are ", num_control_samples, " control samples in the input tables, however, only 1 is expected."))
    }

    write_tsv(results, paste0(lab_batch, "_all_results.tsv"))

    write_tsv(results %>% filter(is_control_sample), paste0(lab_batch, "_control_results.tsv"))

    results_pivoted <- results %>% select(-lab_batch, -is_control_sample) %>% pivot_longer(!sample_id, names_to=c("condition",".value"), names_pattern="([^_]+)_(.+)")
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

    write_tsv(results_summarised, paste0(lab_batch, "_summarised_results.tsv"))

    ggplot(results_pivoted, aes(x=adjusted)) +
      geom_density(aes(color=condition), fill=NA, position = "identity") +
      xlim(-5,5) + theme_bw() + xlab("z-score") + geom_function(fun=dnorm) +
      ylab("density")
    ggsave(filename = paste0(lab_batch, "_score_distribution.png"), dpi=300, width = 6, height = 6)

    writeLines(lab_batch, "lab_batch.txt")

    missing_sites_shifts <-  c("~{sep='","' missing_sites_shifts}") %>% map(read_tsv) %>% reduce(bind_rows)
    write_tsv(missing_sites_shifts, paste0(lab_batch, "_missing_sites_shifts.tsv"))

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File batch_all_results = "~{lab_batch}_all_results.tsv"
    File batch_control_results = "~{lab_batch}_control_results.tsv"
    File batch_summarised_results = "~{lab_batch}_summarised_results.tsv"
    File batch_score_distribution = "~{lab_batch}_score_distribution.png"
    File batch_missing_sites_shifts = "~{lab_batch}_missing_sites_shifts.tsv"
  }
}

task PlotPCA {
  input {
    Array[File] target_pc_projections
    File population_pc_projections
    String lab_batch
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
      geom_point(data=target_pcs, aes(color="~{lab_batch}")) +
      theme_bw()

    ggsave(filename = "~{lab_batch}_PCA_plot.png", dpi=300, width = 6, height = 6)

    EOF

  >>>

  runtime {
    docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File pc_plot = "~{lab_batch}_PCA_plot.png"
  }
}

task BuildHTMLReport {
  input {
    File batch_all_results
    File batch_control_results
    File batch_missing_sites_shifts
    File expected_control_results
    File batch_summarised_results
    File score_distribution
    Array[File] target_pc_projections
    File population_pc_projections
    String population_name
    String lab_batch
  }

  command <<<
    set -xeo pipefail

    cat << EOF > ~{lab_batch}_report.Rmd
    ---
    title: "Batch ~{lab_batch} PRS Summary"
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
    library(tibble)
    library(plotly)

    batch_all_results <- read_tsv("~{batch_all_results}")
    batch_control_results <- read_tsv("~{batch_control_results}")
    expected_control_results <- read_csv("~{expected_control_results}")
    batch_summary <- read_tsv("~{batch_summarised_results}")
    batch_summary <- batch_summary %>% rename_with(.cols = -condition, ~ str_to_title(gsub("_"," ", .x)))
    \`\`\`


    ## Control Sample
    \`\`\`{r control, echo = FALSE, results = "asis"}
    kable(list(batch_control_results, expected_control_results) %>% reduce(bind_rows) %>% select(ends_with('_adjusted')) %>% add_column(sample=c('batch_control', 'expected_control'), .before=1), digits = 2)
    \`\`\`

    ## Batch Summary
    \`\`\`{r summary table, echo = FALSE, results = "asis" }
    kable(batch_summary, digits = 2)
    \`\`\`



    ## Batch Score distribution
    ![](~{score_distribution})

    ## PCA
    #### Hover for sample ID
    \`\`\`{r pca plot, echo=FALSE, message=FALSE, warning=FALSE, results="asis"}
    target_pcs <- c("~{sep='","' target_pc_projections}") %>% map(read_tsv) %>% reduce(bind_rows)
    population_pcs <- read_tsv("~{population_pc_projections}")

    p <- ggplot(population_pcs, aes(x=PC1, y=PC2, color="~{population_name}")) +
      geom_point(size=0.1, alpha=0.1) +
      geom_point(data=target_pcs, aes(color="~{lab_batch}", text=paste0("Sample ID: ", IID))) +
      theme_bw()
    ggplotly(p, tooltip="text")
    \`\`\`

    ## Individual Sample Results (without control sample)
    \`\`\`{r sample results , echo = FALSE, results = "asis"}
    kable(batch_all_results %>% filter(!is_control_sample) %>% select(-is_control_sample) %>% mutate(across(ends_with("risk"), ~ kableExtra::cell_spec(.x, color=ifelse(is.na(.x), "blue", ifelse(.x=="NOT_RESULTED", "red", ifelse(.x == "HIGH", "orange", "green")))))), escape = FALSE, digits = 2, format = "pandoc")
    \`\`\`

    ## Missing sites
    \`\`\`{r missing sites load, include = FALSE}
    batch_missing_sites <- read_tsv("~{batch_missing_sites_shifts}")
    batch_missing_sites <- batch_missing_sites %>% filter(n_missing_sites > 0)
    \`\`\`
    \`r if (batch_missing_sites %>% count() == 0) {"All expected sites were included in all scores for all samples."} else {"Scores missing expected sites are shown below."}\`

    \`\`\`{r missing sites table, echo = FALSE, results = "asis"}
    if (batch_missing_sites %>% count() > 0) {
    kable(batch_missing_sites %>% mutate(across(all_of(c("sample_id", "condition")), ~kableExtra::cell_spec(.x, color = ifelse(pmax(abs(potential_high_percentile - percentile), abs(potential_low_percentile - percentile)) > 0.05, "red", ifelse(pmax(abs(potential_high_percentile - percentile), abs(potential_low_percentile - percentile)) > 0.02, "orange", "black"))))), escape = FALSE, digits = 2, format = "pandoc")    }
    \`\`\`
    EOF

    Rscript -e "library(rmarkdown); rmarkdown::render('~{lab_batch}_report.Rmd', 'html_document')"
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/tidyverse_kableextra_docker@sha256:fd21f5608a3d43add02f8a8490e49db67f078cb2b906f8cd959a9767350b8c24"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File report = "~{lab_batch}_report.html"
  }
}
