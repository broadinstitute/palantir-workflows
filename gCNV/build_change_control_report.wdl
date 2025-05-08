version 1.0

workflow build_change_control_report {
    input {
        Array[File] ppv_files
        Array[File] recall_files
        Array[File] curated_events_match_files
        Array[String] versions
        Array[String] descriptions
        String report_basename
    }

    call build_report {
        input:
            ppv_files = ppv_files,
            recall_files = recall_files,
            curated_events_match_files = curated_events_match_files,
            versions = versions,
            descriptions = descriptions,
            report_basename = report_basename
    }

    output {
        File report = build_report.report
    } 
}

task build_report {
    input {
        Array[File] ppv_files
        Array[File] recall_files
        Array[File] curated_events_match_files
        Array[String] versions
        Array[String] descriptions
        String report_basename
        Int mem_gb = 4
    }

    command <<<
        set -eo pipefail
        Rscript -e "install.packages('kableExtra')"
        tlmgr install booktabs multirow wrapfig float colortbl pdflscape tabu threeparttable threeparttablex ulem makecell xcolor bookmark environ
        cat << EOF > ~{report_basename}_change_control_report.Rmd
        ---
        title: "BGE CNV Change Control"
        date: "\`r Sys.Date()\`"
        output: 
            pdf_document:
                keep_tex: true
        ---

        \`\`\`{r setup, include=FALSE}
        library(dplyr)
        library(ggplot2)
        library(knitr)
        library(readr)
        library(kableExtra)
        library(tidyr)
        library(purrr)
        knitr::opts_chunk\$set(echo = FALSE)
        versions <- tibble(Version=c("~{sep='","' versions}"), Description=c("~{sep='","' descriptions}"))

        ppv_vcfs <- c("~{sep='","' ppv_files}")
        ppv <- ppv_vcfs %>% map(read_tsv) %>% reduce(bind_rows) %>% rowwise() %>% mutate(ppv_low_95=binom.test(round(TP_frac),round(TP_frac+FP_frac))\$conf.int[[1]],
            ppv_high_95=binom.test(round(TP_frac),round(TP_frac+FP_frac))\$conf.int[[2]])

        recall_vcfs <- c("~{sep='","' recall_files}")
        recall <- recall_vcfs %>% map(read_tsv) %>% reduce(bind_rows) %>% rowwise() %>% mutate(recall_low_95=binom.test(round(TP_frac),round(TP_frac+FN_frac))\$conf.int[[1]],
            recall_high_95=binom.test(round(TP_frac),round(TP_frac+FN_frac))\$conf.int[[2]])

        ppv_total_3_exons <- ppv %>% filter(min_exon_count==3) %>% group_by(svtype_bge,version) %>% summarise(TP_frac=sum(TP_frac),
                                                                    FP_frac=sum(FP_frac)) %>%
        mutate(ppv=TP_frac/(TP_frac+FP_frac)) %>% rowwise() %>% mutate(ppv_low_95=binom.test(round(TP_frac),round(TP_frac+FP_frac))\$conf.int[[1]],
            ppv_high_95=binom.test(round(TP_frac),round(TP_frac+FP_frac))\$conf.int[[2]])

        recall_total_3_exons <- recall %>% filter(min_exon_count==3) %>% group_by(svtype_truth,version) %>% summarise(TP_frac=sum(TP_frac), FN_frac=sum(FN_frac)) %>%
        mutate(recall=TP_frac/(TP_frac+FN_frac)) %>% rowwise() %>% mutate(recall_low_95=binom.test(round(TP_frac),round(TP_frac+FN_frac))\$conf.int[[1]],
            recall_high_95=binom.test(round(TP_frac),round(TP_frac+FN_frac))\$conf.int[[2]])
          
        curated_events_files <- c("~{sep='","' curated_events_match_files}")
        curated_events <- curated_events_files %>% map(read_tsv) %>% reduce(bind_rows) %>% select(-overlapping_bge_len)

        found_counts <- curated_events %>% group_by(version) %>% 
        summarise(n_found = sum(bge_overlap_frac > 0.5), n=n()) %>% mutate(found_str = paste0(n_found,"/",n))
        
        curated_events_wide <- curated_events %>% pivot_wider(names_from = version, values_from = bge_overlap_frac, names_prefix = "Observed Frac ") %>% rename(\`Exonic Length\` = curated_event_exon_length) %>% mutate(Event=paste0(contig,":",start,":",\`CNV Type\`)) %>%
        select(-c(\`Patient ID\`, contig, start,\`CNV Type\`)) %>% select(Event, everything())
        \`\`\`

        \`\`\`{r versions}
        kable(versions)
        \`\`\`

        ## PPV and Recall Analysis Description

        This section compares PPV and Recall of the two versions.  GATK-SV calls on 30x WGS data are taken as truth.

        For each event in the BGE validation data, a \`True Positive Fraction\` is awarded corresponding to the fraction of the exon region of the event covered by passing events of the same type in the corresponding truth sample.  A \`False Positive Fraction\` is awarded corresponding to the fraction of the exon region of the event that is not covered by any event (passing or not) of the same type in the corresponding truth sample.  Note that the \`True Positive Fraction\` and \`False Positive Fraction\` may not sum to one if an event is covered by non-passing events in the truth.  In this situation, there is some amount of evidence for the event in the truth, but not enough to make a confident call, so the corresponding portion of the event is considered neither \`True Positive\` nor \`False Positive\`.  The sum of all \`True Positive Fractions\` (\`TP_total\`) and \`False Positive Fractions\` (\`FP_total\`) are calculated, and \`PPV\` is computed as \`TP_total/(TP_total + FP_total)\`.

        For each event in the truth data, a \`True Positive Fraction\` is awarded corresponding to the fraction of the exon region of the event covered by passing events of the same type in the corresponding BGE sample.  A \`False Negative Fraction\` is awarded corresponding to the fraction of the exon region of the event that is not covered by passing events of the same type in the corresponding BGE sample.  Note that for each event, the \`True Positive Fraction\` and \`False Negative Fraction\` will sum to one.  The sum of all \`True Positive Fractions\` (\`TP_total\`) and \`False Negative Fractions\` (\`FP_total\`) are calculated over only events with allele frequency below 1%, \`Recall\` is computed as \`TP_total/(TP_total + FN_total)\`.

        ## PPV as a function of Minimum Number of Exons

        \`\`\`{r ppv_vs_min_exons}
        ggplot(ppv, aes(x=min_exon_count, y=ppv, color=version)) +
        geom_ribbon(aes(ymin=ppv_low_95, ymax=ppv_high_95, fill=version), alpha=0.2, color=NA) +
        geom_line() + 
        facet_grid(svtype_bge~input_material_type_bge) + 
        ylim(0,1) + theme_bw() + xlab("Minimum Number of Exons") +
        scale_x_continuous(breaks = seq(0,10,2)) 
        \`\`\`

        ## Recall as a function of Minimum Number of Exons

        \`\`\`{r recall_vs_min_exons}
        ggplot(recall, aes(x=min_exon_count, y=recall, color=version)) +
        geom_ribbon(aes(ymin=recall_low_95, ymax=recall_high_95, fill=version), alpha=0.2, color=NA) +
        geom_line() + 
        facet_grid(svtype_truth~input_material_type_truth) + 
        ylim(0,1) + theme_bw() + xlab("Minimum Number of Exons") +
        scale_x_continuous(breaks = seq(0,10,2)) 
        \`\`\`


        ## PPV and Recall for 3+ Exons
        \`\`\`{r}
        ppv_recall_3_exons <- ppv_total_3_exons %>% inner_join(recall_total_3_exons, by=c("svtype_bge"="svtype_truth","version")) %>% transmute(svtype_bge, version, ppv_str=paste0(format(ppv, digits=2)," \\\\scriptsize (",format(ppv_low_95, digits=2),
                                                                            "-",format(ppv_high_95,digits=2),")"),
                                        recall_str=paste0(format(recall, digits=2)," \\\\scriptsize (",format(recall_low_95, digits=2),
                                                                            "-",format(recall_high_95,digits=2),")")
                )

        kable(ppv_recall_3_exons%>% mutate(version = gsub("_","\\\\\\\\_",version)), format = "latex", escape=FALSE, booktabs=TRUE, col.names = c("","","ppv \\\\scriptsize (95\\\\% CI)", "recall \\\\scriptsize (95\\\\% CI)"))  %>% collapse_rows(columns = 1)
        \`\`\`

        ## Analytical Sensitivity for Previously Identified Variants in Clinical Samples
        For previously identified variants in clinical samples, \`Observed Frac\` was calculated identically to how fractional \`True Positives\` were calculated for PPV and Recall calculations.  That is, \`Observed Frac\` is the fraction of the exonic regions of the event which is covered by events of the same type in the corresponding BGE sample.  An event was counted as "found" in the BGE data if \`Observed Frac\` was greater than 0.5 for the event.

        \`\`\`{r}
        kable(found_counts %>% select(-c(n_found,n)), col.names = c("", "Events Found in BGE")) 
        \`\`\`

        \`\`\`{r}
        kable(curated_events_wide, booktabs=TRUE, digits = 2) %>%
        column_spec(2, width = "10em") %>% kable_styling(latex_options = c("striped")) %>%
        column_spec(c(4,5), width="7em")
        \`\`\`
        EOF

        Rscript -e "library(rmarkdown); rmarkdown::render('~{report_basename}_change_control_report.Rmd', 'pdf_document')"

    >>>

    runtime {
        docker: "rocker/verse:4.4"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
      }

    output {
        File report_tex = "~{report_basename}_change_control_report.tex"
        File report = "~{report_basename}_change_control_report.pdf"
    }
}