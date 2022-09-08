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
      batch_control_results = AggregateResults.batch_control_results,
      batch_missing_sites_shifts = AggregateResults.batch_missing_sites_shifts,
      expected_control_results = expected_control_results,
      batch_summarised_results = AggregateResults.batch_summarised_results,
      batch_pivoted_results = AggregateResults.batch_pivoted_results,
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

    results_pivoted <- results %>% filter(!is_control_sample) %>% pivot_longer(!c(sample_id, lab_batch, is_control_sample), names_to=c("condition",".value"), names_pattern="([^_]+)_(.+)")
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

    write_tsv(results_pivoted, paste0(lab_batch, "_pivoted_results.tsv"))

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
    File batch_pivoted_results = "~{lab_batch}_pivoted_results.tsv"
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
    File batch_control_results
    File batch_missing_sites_shifts
    File expected_control_results
    File batch_summarised_results
    File batch_pivoted_results
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
    library(DT)
    library(stringi)
    library(tidyr)

    batch_control_results <- read_tsv("~{batch_control_results}", col_types = cols(.default = 'n'))
    expected_control_results <- read_csv("~{expected_control_results}", col_types = cols(.default = 'n'))
    batch_pivoted_results <- read_tsv("~{batch_pivoted_results}")
    batch_summary <- read_tsv("~{batch_summarised_results}")
    batch_summary <- batch_summary %>% rename_with(.cols = -condition, ~ str_to_title(gsub("_"," ", .x)))
    multi_high_samples <- batch_pivoted_results %>% filter(risk=="HIGH") %>% group_by(sample_id) %>%
      summarise(\`high risk conditions\` = paste(condition, collapse = ","), n=n()) %>%
      filter(n>1) %>% select(-n)
    \`\`\`

    \`\`\`{css, echo=FALSE}
    .main-container {
    max-width: 100%;
    margin: auto;
    }

    .plotly {
    margin: auto;
    }
    \`\`\`


    ## Control Sample
    \`\`\`{r control, echo = FALSE, results = "asis", warning = FALSE}
    control_and_expected <- bind_rows(list(batch_control_results, expected_control_results)) %>% select(ends_with('_adjusted'))
    delta_frame_colored <- (control_and_expected[-1,] - control_and_expected[-nrow(control_and_expected),]) %>% mutate(across(everything(), ~ round(.x, digits=2))) %>% mutate(across(everything(), ~ kableExtra::cell_spec(.x, color=ifelse(is.na(.x) || abs(.x) > 0.12, "red", "green"))))
    control_and_expected_char <- control_and_expected %>% mutate(across(everything(), ~ format(round(.x, digits=2), nsmall=2)))
    control_table <- bind_rows(list(control_and_expected_char, delta_frame_colored)) %>% select(order(colnames(.)))
    kable(control_table %>% add_column(sample=c('batch_control', 'expected_control', 'delta'), .before=1), escape = FALSE, digits = 2, format = "pandoc")
    \`\`\`

    ## Batch Summary
    \`\`\`{r summary table, echo = FALSE, results = "asis" }
    kable(batch_summary, digits = 2, escape = FALSE, format = "pandoc")
    \`\`\`

    ## Samples High Risk for Multiple Conditions
    \`r if (multi_high_samples %>% count() == 0) {"No Samples were high risk for multiple conditions."} else {"The following samples were high risk for multiple conditions ."}\`
    \`\`\`{r multi high samples table, echo = FALSE, results = "asis" }
    if (multi_high_samples %>% count() > 0) {
    kable(multi_high_samples, digits = 2, escape = FALSE, format = "pandoc") }
    \`\`\`

    ## Batch Score distribution
    \`\`\`{r score distributions, echo=FALSE, message=FALSE, warning=FALSE, results="asis", fig.align='center'}
    conditions_with_more_than_3_samples <- batch_pivoted_results %>% group_by(condition) %>% filter(!is.na(adjusted)) %>% count() %>% filter(n>3) %>% pull(condition)
    ggplot(batch_pivoted_results %>% filter(condition %in% conditions_with_more_than_3_samples), aes(x=adjusted)) +
      geom_density(aes(color=condition), fill=NA, position = "identity") +
      xlim(-5,5) + theme_bw() + xlab("z-score") + geom_function(fun=dnorm) +
      geom_point(data = batch_pivoted_results %>% filter(!(condition %in% conditions_with_more_than_3_samples)), aes(color=condition, x = adjusted), y=0) +
      ylab("density")
    \`\`\`

    ## PCA
    #### Hover for sample ID
    \`\`\`{r pca plot, echo=FALSE, message=FALSE, warning=FALSE, results="asis", fig.align='center'}
    target_pcs <- c("~{sep='","' target_pc_projections}") %>% map(read_tsv) %>% reduce(bind_rows)
    population_pcs <- read_tsv("~{population_pc_projections}")

    p <- ggplot(population_pcs, aes(x=PC1, y=PC2, color="~{population_name}")) +
      geom_point() +
      geom_point(data=target_pcs, aes(color="~{lab_batch}", text=paste0("Sample ID: ", IID))) +
      theme_bw()
    ggplotly(p, tooltip="text")
    \`\`\`

    ## Individual Sample Results (without control sample)
    \`\`\`{r sample results , echo = FALSE, results = "asis", warning = FALSE, message = FALSE}
    batch_high_counts_per_sample <- batch_pivoted_results %>% group_by(sample_id) %>% summarise(n_high_risk = sum(ifelse(!is.na(risk) & risk =="HIGH", 1, 0)))
    batch_results_table <- batch_pivoted_results %>% filter(!is_control_sample) %>% select(!is_control_sample) %>%
      mutate(across(!c(sample_id, lab_batch, reason_not_resulted, condition), ~kableExtra::cell_spec(gsub("_", " ", ifelse(is.na(as.numeric(.x)), ifelse(is.na(.x), 'SCORE NOT REQUESTED', .x), round(as.numeric(.x), 2))), color=ifelse(is.na(risk), "lightgrey", ifelse(risk=="NOT_RESULTED", "red", ifelse(risk == "HIGH", "orange", "green")))))) %>% # round numbers, color all by risk
      mutate(reason_not_resulted = ifelse(is.na(reason_not_resulted), reason_not_resulted, kableExtra::cell_spec(reason_not_resulted, color="red"))) %>% # reason not resulted always red if exists
      pivot_wider(id_cols = c(sample_id, lab_batch), names_from = condition, names_glue = "{condition}_{.value}", values_from = c(raw, adjusted, percentile, risk, reason_not_resulted)) %>% # pivot to wide format
      inner_join(batch_high_counts_per_sample) # add number of high risk conditions for each sample

    #order columns as desired
    cols <- batch_results_table %>% select(-sample_id, -lab_batch, -n_high_risk) %>% colnames()
    desired_order_values <- c("raw", "adjusted", "percentile", "risk", "reason_not_resulted")
    col_order <- c("sample_id", "lab_batch", "n_high_risk", cols[order(sapply(stri_split_fixed(cols, "_", n=2), "[",1), match(sapply(stri_split_fixed(cols, "_", n=2), "[",2), desired_order_values))])
    batch_results_table <- batch_results_table %>% select(all_of(col_order)) %>%
      rename_with(.cols = ends_with("percentile"), .fn = ~gsub("_percentile", " %", .x,fixed=TRUE)) %>%
      rename_with(.cols = ends_with("adjusted"), .fn = ~gsub("_adjusted", "_adj", .x,fixed=TRUE))

    all_cols = batch_results_table %>% colnames()
    risk_cols = which(endsWith(all_cols, "risk"))
    raw_cols = which(endsWith(all_cols, "raw"))
    adjusted_cols = which(endsWith(all_cols, "adj"))
    percentile_cols = which(endsWith(all_cols, "%"))
    reason_not_resulted_cols = which(endsWith(all_cols, "reason_not_resulted"))
    numeric_cols = batch_results_table %>% select(where(is.numeric)) %>% colnames()
    DT::datatable(batch_results_table,
                  escape = FALSE,
                  extensions = c("Buttons", "FixedColumns"),
                  rownames = FALSE,
                  options=list(scrollX = TRUE,
                  fixedColumns = list(leftColumns = 2),
                  pageLength = 20,
                  lengthMenu = c(10, 20, 50, 100),
                  dom = 'Blfrtip',
                  buttons = list(
                                list(
                                  extend = 'columnToggle',
                                  columns = raw_cols - 1,
                                  text = "Raw Scores"
                                ),
                                list(
                                  extend = 'columnToggle',
                                  columns = adjusted_cols - 1,
                                  text = "Adjusted Scores"
                                ),
                                list(
                                  extend = 'columnToggle',
                                  columns = percentile_cols - 1,
                                  text = "Percentile"
                                ),
                                list(
                                  extend = 'columnToggle',
                                  columns = risk_cols - 1,
                                  text = "Risk"
                                ),
                                list(
                                  extend = 'columnToggle',
                                  columns = reason_not_resulted_cols - 1,
                                  text = "Reason Not Resulted"
                                )
                          )
                  )
    )  %>%
    formatStyle(columns = c("sample_id", "lab_batch"), fontWeight = 'bold')
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
    docker: "us.gcr.io/broad-dsde-methods/tidyverse_kableextra_docker@sha256:f9ad840130f45cabe53d2464e3d5fc4130fd8964e263bab9b3de79d45021e1a1"
    disks: "local-disk 100 HDD"
    memory: "4 GB"
  }

  output {
    File report = "~{lab_batch}_report.html"
  }
}
