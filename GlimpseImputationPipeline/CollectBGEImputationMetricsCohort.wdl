version 1.0

workflow CollectBGEImputationMetricsCohort {
    input {
        Array[String] ancestries
        Boolean collect_af_0_1_single_number = false

        Array[String] sample_ids1
        File eval_vcf1
        Array[String] truth_sample_ids1
        File truth_vcf1
        String configuration_label1
        File annotation_vcf1
        Map[String, String] ancestry_to_af_annotation_map1
        String? intervals1

        Array[String]? sample_ids2
        File? eval_vcf2
        Array[String]? truth_sample_ids2
        File? truth_vcf2
        String? configuration_label2
        File? annotation_vcf2
        Map[String, String] ancestry_to_af_annotation_map2
        String? intervals2

        Array[String]? sample_ids3
        File? eval_vcf3
        Array[String]? truth_sample_ids3
        File? truth_vcf3
        String? configuration_label3
        File? annotation_vcf3
        Map[String, String] ancestry_to_af_annotation_map3
        String? intervals3

        Array[String]? sample_ids4
        File? eval_vcf4
        Array[String]? truth_sample_ids4
        File? truth_vcf4
        String? configuration_label4
        File? annotation_vcf4
        Map[String, String] ancestry_to_af_annotation_map4
        String? intervals4

        Array[String]? sample_ids5
        File? eval_vcf5
        Array[String]? truth_sample_ids5
        File? truth_vcf5
        String? configuration_label5
        File? annotation_vcf5
        Map[String, String] ancestry_to_af_annotation_map5
        String? intervals5

        Array[String]? sample_ids6
        File? eval_vcf6
        Array[String]? truth_sample_ids6
        File? truth_vcf6
        String? configuration_label6
        File? annotation_vcf6
        Map[String, String] ancestry_to_af_annotation_map6
        String? intervals6

        Array[String]? sample_ids7
        File? eval_vcf7
        Array[String]? truth_sample_ids7
        File? truth_vcf7
        String? configuration_label7
        File? annotation_vcf7
        Map[String, String] ancestry_to_af_annotation_map7
        String? intervals7

        String output_basename = "cohort"
        Int plot_width = 14
        Int plot_height = 5

        Int preemptible = 1
    }

    call PearsonCorrelationByAF as PearsonByAF {
        input:
            evalVcf = eval_vcf1,
            af_resource = annotation_vcf1,
            ancestries = ancestries,
            ancestry_to_af_annotation_map = ancestry_to_af_annotation_map1,
            truthVcf = truth_vcf1,
            intervals = intervals1,
            output_basename = output_basename,
            eval_sample_ids = sample_ids1,
            truth_sample_ids = truth_sample_ids1,
            preemptible = preemptible
    }

    if (collect_af_0_1_single_number) {
        call PearsonCorrelationByAF as PearsonByAF01 {
            input:
                evalVcf = eval_vcf1,
                af_resource = annotation_vcf1,
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map1,
                truthVcf = truth_vcf1,
                min_af_for_accuracy_metrics = 0.1,
                n_bins = 2,
                right_edge_first_bin = 0.1,
                intervals = intervals1,
                output_basename = output_basename,
                eval_sample_ids = sample_ids1,
                truth_sample_ids = truth_sample_ids1,
                preemptible = preemptible
        }
    }

    if (defined(eval_vcf2)) {
        call PearsonCorrelationByAF as PearsonByAF_2 {
            input:
                evalVcf = select_first([eval_vcf2]),
                af_resource = select_first([annotation_vcf2]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map2,
                truthVcf = select_first([truth_vcf2]),
                intervals = intervals2,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids2]),
                truth_sample_ids = select_first([truth_sample_ids2]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_2 {
                input:
                    evalVcf = select_first([eval_vcf2]),
                    af_resource = select_first([annotation_vcf2]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map2,
                    truthVcf = select_first([truth_vcf2]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals2,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids2]),
                    truth_sample_ids = select_first([truth_sample_ids2]),
                    preemptible = preemptible
            }
        }
    }

    if (defined(eval_vcf3)) {
        call PearsonCorrelationByAF as PearsonByAF_3 {
            input:
                evalVcf = select_first([eval_vcf3]),
                af_resource = select_first([annotation_vcf3]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map3,
                truthVcf = select_first([truth_vcf3]),
                intervals = intervals3,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids3]),
                truth_sample_ids = select_first([truth_sample_ids3]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_3 {
                input:
                    evalVcf = select_first([eval_vcf3]),
                    af_resource = select_first([annotation_vcf3]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map3,
                    truthVcf = select_first([truth_vcf3]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals3,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids3]),
                    truth_sample_ids = select_first([truth_sample_ids3]),
                    preemptible = preemptible
            }
        }
    }

    if (defined(eval_vcf4)) {
        call PearsonCorrelationByAF as PearsonByAF_4 {
            input:
                evalVcf = select_first([eval_vcf4]),
                af_resource = select_first([annotation_vcf4]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map4,
                truthVcf = select_first([truth_vcf4]),
                intervals = intervals4,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids4]),
                truth_sample_ids = select_first([truth_sample_ids4]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_4 {
                input:
                    evalVcf = select_first([eval_vcf4]),
                    af_resource = select_first([annotation_vcf4]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map4,
                    truthVcf = select_first([truth_vcf4]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals4,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids4]),
                    truth_sample_ids = select_first([truth_sample_ids4]),
                    preemptible = preemptible
            }
        }
    }

    if (defined(eval_vcf5)) {
        call PearsonCorrelationByAF as PearsonByAF_5 {
            input:
                evalVcf = select_first([eval_vcf5]),
                af_resource = select_first([annotation_vcf5]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map5,
                truthVcf = select_first([truth_vcf5]),
                intervals = intervals5,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids5]),
                truth_sample_ids = select_first([truth_sample_ids5]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_5 {
                input:
                    evalVcf = select_first([eval_vcf3]),
                    af_resource = select_first([annotation_vcf5]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map5,
                    truthVcf = select_first([truth_vcf5]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals5,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids5]),
                    truth_sample_ids = select_first([truth_sample_ids5]),
                    preemptible = preemptible
            }
        }
    }

    if (defined(eval_vcf6)) {
        call PearsonCorrelationByAF as PearsonByAF_6 {
            input:
                evalVcf = select_first([eval_vcf6]),
                af_resource = select_first([annotation_vcf6]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map6,
                truthVcf = select_first([truth_vcf6]),
                intervals = intervals6,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids6]),
                truth_sample_ids = select_first([truth_sample_ids6]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_6 {
                input:
                    evalVcf = select_first([eval_vcf3]),
                    af_resource = select_first([annotation_vcf6]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map6,
                    truthVcf = select_first([truth_vcf6]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals6,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids6]),
                    truth_sample_ids = select_first([truth_sample_ids6]),
                    preemptible = preemptible
            }
        }
    }

    if (defined(eval_vcf7)) {
        call PearsonCorrelationByAF as PearsonByAF_7 {
            input:
                evalVcf = select_first([eval_vcf7]),
                af_resource = select_first([annotation_vcf7]),
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map7,
                truthVcf = select_first([truth_vcf7]),
                intervals = intervals7,
                output_basename = output_basename,
                eval_sample_ids = select_first([sample_ids7]),
                truth_sample_ids = select_first([truth_sample_ids7]),
                preemptible = preemptible
        }

        if (collect_af_0_1_single_number) {
            call PearsonCorrelationByAF as PearsonByAF01_7 {
                input:
                    evalVcf = select_first([eval_vcf3]),
                    af_resource = select_first([annotation_vcf7]),
                    ancestries = ancestries,
                    ancestry_to_af_annotation_map = ancestry_to_af_annotation_map7,
                    truthVcf = select_first([truth_vcf7]),
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals7,
                    output_basename = output_basename,
                    eval_sample_ids = select_first([sample_ids7]),
                    truth_sample_ids = select_first([truth_sample_ids7]),
                    preemptible = preemptible
            }
        }
    }

    call GenerateCorrelationPlots {
        input:
            correlation_file1 = PearsonByAF.correlations,
            correlation_file2 = PearsonByAF_2.correlations,
            correlation_file3 = PearsonByAF_3.correlations,
            correlation_file4 = PearsonByAF_4.correlations,
            correlation_file5 = PearsonByAF_5.correlations,
            correlation_file6 = PearsonByAF_6.correlations,
            correlation_file7 = PearsonByAF_7.correlations,
            ancestries = ancestries,
            eval_sample_ids1 = sample_ids1,
            eval_sample_ids2 = sample_ids2,
            eval_sample_ids3 = sample_ids3,
            eval_sample_ids4 = sample_ids4,
            eval_sample_ids5 = sample_ids5,
            eval_sample_ids6 = sample_ids6,
            eval_sample_ids7 = sample_ids7,
            configuration_label1 = configuration_label1,
            configuration_label2 = configuration_label2,
            configuration_label3 = configuration_label3,
            configuration_label4 = configuration_label4,
            configuration_label5 = configuration_label5,
            configuration_label6 = configuration_label6,
            configuration_label7 = configuration_label7,
            plot_width = plot_width,
            plot_height = plot_height,
            preemptible = preemptible
    }

    output {
        File correlation_plot = GenerateCorrelationPlots.correlation_plot
        File correlation_data = GenerateCorrelationPlots.correlation_data
        File? correlation_file1_0_1 = PearsonByAF01.correlations
        File? correlation_file2_0_1 = PearsonByAF01_2.correlations
        File? correlation_file3_0_1 = PearsonByAF01_3.correlations
        File? correlation_file4_0_1 = PearsonByAF01_4.correlations
        File? correlation_file5_0_1 = PearsonByAF01_5.correlations
        File? correlation_file6_0_1 = PearsonByAF01_6.correlations
        File? correlation_file7_0_1 = PearsonByAF01_7.correlations
    }
}



task GenerateCorrelationPlots {
    input {
        File correlation_file1
        File? correlation_file2
        File? correlation_file3
        File? correlation_file4
        File? correlation_file5
        File? correlation_file6
        File? correlation_file7
        Array[String] ancestries
        Array[String] eval_sample_ids1
        Array[String]? eval_sample_ids2
        Array[String]? eval_sample_ids3
        Array[String]? eval_sample_ids4
        Array[String]? eval_sample_ids5
        Array[String]? eval_sample_ids6
        Array[String]? eval_sample_ids7
        String configuration_label1
        String? configuration_label2
        String? configuration_label3
        String? configuration_label4
        String? configuration_label5
        String? configuration_label6
        String? configuration_label7

        Int plot_width
        Int plot_height

        Int mem_gb = 2
        Int preemptible = 1
    }

    command <<<
        python <<'EOF'

import pandas as pd

correlation_dfs = []

ancestries = ['~{sep="', '" ancestries}']
eval_sample_ids1 = ['~{sep="', '" eval_sample_ids1}']
ancestry_dict = dict()
for i in range(len(eval_sample_ids1)):
    ancestry_dict[eval_sample_ids1[i]] = ancestries[i]

df = pd.read_csv('~{correlation_file1}', sep="\t", comment='#')
if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
    raise RuntimeError(f'~{correlation_file1} has no sites in it.')
df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
df['configuration'] = '~{configuration_label1}'
correlation_dfs.append(df)

if ~{if defined(correlation_file2) then "True" else "False"}:
    eval_sample_ids2 = ['~{sep="', '" eval_sample_ids2}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids2)):
        ancestry_dict[eval_sample_ids2[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file2}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file2} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label2}'
    correlation_dfs.append(df)

if ~{if defined(correlation_file3) then "True" else "False"}:
    eval_sample_ids3 = ['~{sep="', '" eval_sample_ids3}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids3)):
        ancestry_dict[eval_sample_ids3[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file3}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file3} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label3}'
    correlation_dfs.append(df)

if ~{if defined(correlation_file4) then "True" else "False"}:
    eval_sample_ids4 = ['~{sep="', '" eval_sample_ids4}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids4)):
        ancestry_dict[eval_sample_ids4[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file4}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file4} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label4}'
    correlation_dfs.append(df)

if ~{if defined(correlation_file5) then "True" else "False"}:
    eval_sample_ids5 = ['~{sep="', '" eval_sample_ids5}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids5)):
        ancestry_dict[eval_sample_ids5[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file5}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file5} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label5}'
    correlation_dfs.append(df)

if ~{if defined(correlation_file6) then "True" else "False"}:
    eval_sample_ids6 = ['~{sep="', '" eval_sample_ids6}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids6)):
        ancestry_dict[eval_sample_ids6[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file6}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file6} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label6}'
    correlation_dfs.append(df)

if ~{if defined(correlation_file7) then "True" else "False"}:
    eval_sample_ids7 = ['~{sep="', '" eval_sample_ids7}']
    ancestry_dict = dict()
    for i in range(len(eval_sample_ids7)):
        ancestry_dict[eval_sample_ids7[i]] = ancestries[i]

    df = pd.read_csv('~{correlation_file7}', sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        raise RuntimeError(f'~{correlation_file7} has no sites in it.')
    df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
    df['configuration'] = '~{configuration_label7}'
    correlation_dfs.append(df)

pd.concat(correlation_dfs).to_csv('correlation_data.tsv', sep='\t')
EOF
        cat <<'EOF' > script.R
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(magrittr)
library(svglite)

corr <- read_tsv("correlation_data.tsv")

corr_gathered <- gather(corr, key="type", value="correlation", SNP_CORRELATION, INDEL_CORRELATION) %>%
  mutate(variant_type=ifelse(grepl("SNP", type, fixed=TRUE), "snp", "indel")) %>%
  transform(correlation = as.numeric(correlation))
ggplot(corr_gathered %>% filter(!is.na(correlation)), aes(x=BIN_CENTER, y=correlation^2)) +
  geom_point(size=0.2, alpha=1.0, aes(color=configuration)) +
  #geom_line(aes(color=ancestry)) +
  #geom_smooth(se=FALSE, aes(color=branch)) +
  geom_smooth(aes(color=configuration)) +
  facet_grid(variant_type ~ ancestry) + scale_x_log10() + theme_bw() +
  xlab("Minor Allele Frequency") + ylab(bquote(R^2))

ggsave(filename="correlation_plot.svg", width=~{plot_width}, height=~{plot_height})
EOF
        Rscript script.R
    >>>

    output {
        File correlation_plot = "correlation_plot.svg"
        File correlation_data = "correlation_data.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/mgatzen/bge_metrics_collection:latest"
        preemptible: preemptible
        cpu: 2
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
    }

}

task PearsonCorrelationByAF {
    input {
        File evalVcf
        File truthVcf
        Array[String] eval_sample_ids
        Array[String] truth_sample_ids
        Map[String, String] ancestry_to_af_annotation_map
        String output_basename
        Array[String] ancestries
        File af_resource
        File? sites
        String? intervals
        String? dosage_field
        Int? n_bins
        Float? right_edge_first_bin
        Float? min_af_for_accuracy_metrics
        Int mem_gb = 16
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail

        python <<'EOF'
eval_sample_ids = ['~{sep="', '" eval_sample_ids}']
truth_sample_ids = ['~{sep="', '" truth_sample_ids}']
ancestries = ['~{sep="', '" ancestries}']

# Write sample-map (eval to truth mapping)
with open('sample_map.list', 'w') as sample_map_file:
    for i in range(len(eval_sample_ids)):
        sample_map_file.write(f'{eval_sample_ids[i]}:{truth_sample_ids[i]}\n')

# Write af_expressions (truth to ancestry-specific af annotation mapping)
ancestry_to_af_annotation_dict = dict()
with open('~{write_map(ancestry_to_af_annotation_map)}') as map_file:
    for line in map_file:
        ancestry, af_annotation = line.strip().split('\t')
        ancestry_to_af_annotation_dict[ancestry] = af_annotation

with open('af_expressions.list', 'w') as af_expressions_file:
    for i in range(len(truth_sample_ids)):
        af_expressions_file.write(f'{truth_sample_ids[i]}:{ancestry_to_af_annotation_dict[ancestries[i]]}\n')
EOF

        gatk --java-options "-Xmx~{mem_gb - 2}G" EvaluateGenotypingPerformance -eval ~{evalVcf} -truth ~{truthVcf} --af-annotations af_expressions.list --resource ~{af_resource} \
        ~{"--ids " + sites} ~{"-L " + intervals} --sample-map sample_map.list ~{"--dosage-field " + dosage_field} -O ~{output_basename}.correlations.tsv \
        -OA ~{output_basename}.accuracy.tsv ~{"-nbins " + n_bins} ~{"-first-bin-right-edge " + right_edge_first_bin} ~{"--min-af-for-accuracy-metrics " + min_af_for_accuracy_metrics} --allow-differing-ploidies
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:5910defbc137d43145e6cf1f9c3539ae418b6e465250d55465e19e77773445c4"
        disks: "local-disk 500 HDD"
        memory: mem_gb + " GB"
        preemptible: preemptible
    }

    output {
        File correlations = "~{output_basename}.correlations.tsv"
        File accuracy = "~{output_basename}.accuracy.tsv"
    }
}