version 1.0

workflow PlotROC {
    input {
        String dataset_name
        Array[File] files
        Int? preemptible
    }

    call PlotROCTask {
        input:
            dataset_name = dataset_name,
            files = files,
            preemptible = preemptible
    }

    output {
        File result = PlotROCTask.result
        File snp_plot = PlotROCTask.snp_plot
        File snp_plot_zoomed = PlotROCTask.snp_plot_zoomed
        File snp_plot_zoomed_LCR = PlotROCTask.snp_plot_zoomed_LCR
        File indel_plot = PlotROCTask.indel_plot
        File indel_plot_zoomed = PlotROCTask.indel_plot_zoomed
        File indel_plot_zoomed_LCR = PlotROCTask.indel_plot_zoomed_LCR
        File table = PlotROCTask.table
        File best_fscore = PlotROCTask.best_fscore
    }
}

task PlotROCTask {
    input {
        String dataset_name
        Array[File] files
        Int? preemptible
    }

    command <<<
        for file in ~{sep=' ' files}; do
            cp $file .
        done

        cat <<EOF > script.R
library(tidyverse)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

dataset = args[1]
home = "/cromwell_root/"

files = list.files(home, pattern = "*.tsv", full.names = TRUE)
df_all = tibble()

for (file in files){
  matches=str_match(basename(file), pattern = "([[:alnum:]-]+)_([:alnum:]+)_([:upper:]{3})?_?vcfeval_([a-z_]+)_roc.tsv")
  tool=matches[2]
  data=matches[3]
  region=ifelse(is.na(matches[4]), "all", matches[4])
  type=matches[5]
  df = read_tsv(file, skip = 6, col_types = "dddddddd") %>% rename(score = `#score`) %>% mutate(tool = tool, data = data, region=region, type = type, file = file)
  df_all = bind_rows(df_all, df)
}

### SNPs
p_snp = ggplot(df_all %>% filter(type == "snp")) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) + #, size = 0.001 after tool)
  ggtitle(str_c("snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x") + theme(legend.position = "bottom")
ggsave(filename = str_c(dataset, "_snp.png"), plot = p_snp, width = 8.5, height = 9)
if (df_all %>% filter(type == "snp") %>% filter(sensitivity > 0.93) %>% count() > 0) {
  p_snp_zoomed = ggplot(df_all %>% filter(type == "snp") %>% filter(sensitivity > 0.93)) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) +
    ggtitle(str_c("snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x")  + theme(legend.position = "bottom")
} else {
  p_snp_zoomed = ggplot() + annotate("text", x = 0, y = 0, size=5, label = str_c("No data for ", dataset, " SNPs with sensitivity > 0.93")) + theme_void()
}
ggsave(filename = str_c(dataset, "_snp_zoomed.png"), plot = p_snp_zoomed, width = 8.5, height = 9)
if (df_all %>% filter(type == "snp" & region == "LCR") %>% filter(sensitivity > 0.75) %>% count() > 0) {
  p_snp_zoomed_LCR = ggplot(df_all %>% filter(type == "snp" & region == "LCR") %>% filter(sensitivity > 0.75)) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) +
    ggtitle(str_c("snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x")  + theme(legend.position = "bottom")
} else {
  p_snp_zoomed_LCR = ggplot() + annotate("text", x = 0, y = 0, size=5, label = str_c("No data for ", dataset, " LCR SNPs with sensitivity > 0.75")) + theme_void()
}
ggsave(filename = str_c(dataset, "_snp_zoomed_LCR.png"), plot = p_snp_zoomed_LCR, width = 8.5, height = 9)

### INDELs
p_indel= ggplot(df_all %>% filter(type == "non_snp")) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) +
  ggtitle(str_c("non_snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x") + theme(legend.position = "bottom")
ggsave(filename = str_c(dataset, "_indel.png"), plot = p_indel, width = 8.5, height = 9)
if (df_all %>% filter(type == "non_snp") %>% filter(sensitivity > 0.8) %>% count() > 0) {
  p_indel_zoomed = ggplot(df_all %>% filter(type == "non_snp") %>% filter(sensitivity > 0.8)) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) +
    ggtitle(str_c("non_snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x") + theme(legend.position = "bottom")
} else {
  p_indel_zoomed = ggplot() + annotate("text", x = 0, y = 0, size=5, label = str_c("No data for ", dataset, " INDELs with sensitivity > 0.8")) + theme_void()
}
ggsave(filename = str_c(dataset, "_indel_zoomed.png"), plot = p_indel_zoomed, width = 8.5, height = 9)
if (df_all %>% filter(type == "non_snp" & region == "LCR") %>% filter(sensitivity > 0.7) %>% count() > 0) {
  p_indel_zoomed_LCR = ggplot(df_all %>% filter(type == "non_snp" & region == "LCR") %>% filter(sensitivity > 0.7)) + geom_point(aes(x = false_positives, y = sensitivity, color = tool)) +
    ggtitle(str_c("non_snp: ", dataset)) + facet_grid(. ~ region, scale = "free_x") + theme(legend.position = "bottom")
} else {
  p_indel_zoomed_LCR = ggplot() + annotate("text", x = 0, y = 0, size=5, label = str_c("No data for ", dataset, " LCR INDELs with sensitivity > 0.7")) + theme_void()
}
ggsave(filename = str_c(dataset, "_indel_zoomed_LCR.png"), plot = p_indel_zoomed_LCR, width = 8.5, height = 9)

write_tsv(df_all, str_c(dataset, ".tsv"))

# There should be (2 * 3 * 2 = 12 rows).
best_fscore = df_all %>% group_by(tool, region, type) %>% filter(f_measure == max(f_measure)) %>% 
  filter(true_positives_call == max(true_positives_call)) %>% filter(false_positives == min(false_positives)) %>% mutate(dataset = dataset)
write_tsv(best_fscore, str_c(dataset, "_best_fscore.tsv"))

# Assuming this script is called by a single sample, output the best fscore for SNPs and INDEL to a file.
snp_qual_threshold = best_fscore %>% filter(region == "all" & type == "snp") %>% ungroup() %>% top_n(n = 1, wt = score) %>% select(score) %>% pull()
indel_qual_threshold = best_fscore %>% filter(region == "all" & type == "non_snp") %>% ungroup() %>% top_n(n = 1, wt = score) %>% select(score) %>% pull()

write(snp_qual_threshold, "snp_qual_threshold.txt")
write(indel_qual_threshold, "indel_qual_threshold.txt")
EOF

        Rscript script.R ~{dataset_name}

        ls >> result.txt
        echo "end ls" >> result.txt
    >>>

    output {
        File result = "result.txt"
        File snp_plot = dataset_name + "_snp.png"
        File snp_plot_zoomed = dataset_name + "_snp_zoomed.png"
        File snp_plot_zoomed_LCR = dataset_name + "_snp_zoomed_LCR.png"
        File indel_plot = dataset_name + "_indel.png"
        File indel_plot_zoomed = dataset_name + "_indel_zoomed.png"
        File indel_plot_zoomed_LCR = dataset_name + "_indel_zoomed_LCR.png"
        File table = dataset_name + ".tsv"
        File best_fscore = dataset_name + "_best_fscore.tsv"
    }

    runtime {
        docker: "rocker/tidyverse"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk 200 HDD"
    }
}