version 1.0

workflow CollectBGEImputationMetrics {
    input {
        Array[String] ancestries
        String? intervals

        Array[String] sample_ids1
		Array[File] eval_vcfs1
        Array[String] truth_sample_ids1
		Array[File] truth_vcfs1
        String configuration_label1
		File annotation_vcf1
        Map[String, String] ancestry_to_af_annotation_map1

        Array[String]? sample_ids2
		Array[File]? eval_vcfs2
        Array[String]? truth_sample_ids2
		Array[File]? truth_vcfs2
        String? configuration_label2
		File? annotation_vcf2
        Map[String, String] ancestry_to_af_annotation_map2

        Array[String]? sample_ids3
		Array[File]? eval_vcfs3
        Array[String]? truth_sample_ids3
		Array[File]? truth_vcfs3
        String? configuration_label3
		File? annotation_vcf3
        Map[String, String] ancestry_to_af_annotation_map3
	}

    scatter(ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id in zip(ancestries, zip(zip(truth_vcfs1, truth_sample_ids1), zip(eval_vcfs1, sample_ids1)))) {
        String ancestry1 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.left
        File eval_vcf1 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.left
        String eval_sample_id1 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.right
        File truth_vcf1 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.left
        String truth_sample_id1 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.right

        call PearsonCorrelationByAF as PearsonByAF {
            input:
                evalVcf = eval_vcf1,
                af_resource = annotation_vcf1,
                af_expression = ancestry_to_af_annotation_map1[ancestry1],
                truthVcf = truth_vcf1,
                intervals = intervals,
                output_basename = eval_sample_id1,
                eval_sample_id = eval_sample_id1,
                truth_sample_id = truth_sample_id1,
        }

        call PearsonCorrelationByAF as PearsonByAF01 {
            input:
                evalVcf = eval_vcf1,
                af_resource = annotation_vcf1,
                af_expression = ancestry_to_af_annotation_map1[ancestry1],
                truthVcf = truth_vcf1,
                min_af_for_accuracy_metrics = 0.1,
                n_bins = 2,
                right_edge_first_bin = 0.1,
                intervals = intervals,
                output_basename = eval_sample_id1,
                eval_sample_id = eval_sample_id1,
                truth_sample_id = truth_sample_id1
        }
    }

    if (defined(eval_vcfs2)) {
        scatter(ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id in zip(ancestries, zip(zip(select_first([truth_vcfs2, []]), select_first([truth_sample_ids2, []])), zip(select_first([eval_vcfs2, []]), select_first([sample_ids2, []]))))) {
            String ancestry2 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.left
            File eval_vcf2 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.left
            String eval_sample_id2 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.right
            File truth_vcf2 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.left
            String truth_sample_id2 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.right

            call PearsonCorrelationByAF as PearsonByAF_2 {
                input:
                    evalVcf = eval_vcf2,
                    af_resource = select_first([annotation_vcf2, []]),
                    af_expression = ancestry_to_af_annotation_map2[ancestry2],
                    truthVcf = truth_vcf2,
                    intervals = intervals,
                    output_basename = eval_sample_id2,
                    eval_sample_id = eval_sample_id2,
                    truth_sample_id = truth_sample_id2,
            }

            call PearsonCorrelationByAF as PearsonByAF01_2 {
                input:
                    evalVcf = eval_vcf2,
                    af_resource = select_first([annotation_vcf2, []]),
                    af_expression = ancestry_to_af_annotation_map2[ancestry2],
                    truthVcf = truth_vcf2,
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals,
                    output_basename = eval_sample_id2,
                    eval_sample_id = eval_sample_id2,
                    truth_sample_id = truth_sample_id2
            }
        }
    }

    if (defined(eval_vcfs3)) {
        scatter(ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id in zip(ancestries, zip(zip(select_first([truth_vcfs3, []]), select_first([truth_sample_ids3, []])), zip(select_first([eval_vcfs3, []]), select_first([sample_ids3, []]))))) {
            String ancestry3 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.left
            File eval_vcf3 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.left
            String eval_sample_id3 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.right.right
            File truth_vcf3 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.left
            String truth_sample_id3 = ancestry___and___truth_vcf_and_truth_sample_id__and__eval_vcf_and_eval_sample_id.right.left.right

            call PearsonCorrelationByAF as PearsonByAF_3 {
                input:
                    evalVcf = eval_vcf3,
                    af_resource = select_first([annotation_vcf3, []]),
                    af_expression = ancestry_to_af_annotation_map3[ancestry3],
                    truthVcf = truth_vcf3,
                    intervals = intervals,
                    output_basename = eval_sample_id3,
                    eval_sample_id = eval_sample_id3,
                    truth_sample_id = truth_sample_id3,
            }

            call PearsonCorrelationByAF as PearsonByAF01_3 {
                input:
                    evalVcf = eval_vcf3,
                    af_resource = select_first([annotation_vcf3, []]),
                    af_expression = ancestry_to_af_annotation_map3[ancestry3],
                    truthVcf = truth_vcf3,
                    min_af_for_accuracy_metrics = 0.1,
                    n_bins = 2,
                    right_edge_first_bin = 0.1,
                    intervals = intervals,
                    output_basename = eval_sample_id3,
                    eval_sample_id = eval_sample_id3,
                    truth_sample_id = truth_sample_id3
            }
        }
    }

    call GenerateCorrelationPlots {
        input:
            correlation_files1 = PearsonByAF.correlations,
            correlation_files2 = PearsonByAF_2.correlations,
            correlation_files3 = PearsonByAF_3.correlations,
            ancestries = ancestries,
            configuration_label1 = configuration_label1,
            configuration_label2 = configuration_label2,
            configuration_label3 = configuration_label3
    }

    output {
        File correlation_plot = GenerateCorrelationPlots.correlation_plot
        Array[File] correlation_files1_0_1 = PearsonByAF01.correlations
        Array[File]? correlation_files2_0_1 = PearsonByAF01_2.correlations
        Array[File]? correlation_files3_0_1 = PearsonByAF01_3.correlations
    }
}



task GenerateCorrelationPlots {
    input {
        Array[File] correlation_files1
        Array[File]? correlation_files2
        Array[File]? correlation_files3
        Array[String] ancestries
        String configuration_label1
        String? configuration_label2
        String? configuration_label3

        Int mem_gb = 2
        Int preemptible = 1
    }

    command <<<
        python <<'EOF'

import pandas as pd

correlation_filenames = ['~{sep="', '" correlation_files1}']
ancestries = ['~{sep="', '" ancestries}']

correlation_dfs = []
for i in range(len(correlation_filenames)):
    filename = correlation_filenames[i]
    ancestry = ancestries[i]
    if ancestry not in {'AFR', 'EAS', 'NFE', 'SAS'}:
        continue
    df = pd.read_csv(filename, sep="\t", comment='#')
    if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
        print(f'{filename} has no sites in it.')
        continue
    df['ancestry'] = ancestry
    df['configuration'] = '~{configuration_label1}'
    correlation_dfs.append(df)

correlation_data = pd.concat(correlation_dfs)

if ~{if defined(correlation_files2) then "True" else "False"}:
    correlation_filenames2 = ['~{sep="', '" correlation_files2}']
    for i in range(len(correlation_filenames2)):
        filename = correlation_filenames2[i]
        ancestry = ancestries[i]
        if ancestry not in {'AFR', 'EAS', 'NFE', 'SAS'}:
            continue
        df = pd.read_csv(filename, sep="\t", comment='#')
        if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
            print(f'{filename} has no sites in it.')
            continue
        df['ancestry'] = ancestry
        df['configuration'] = '~{configuration_label2}'
        correlation_dfs.append(df)

    correlation_data2 = pd.concat(correlation_dfs)
    correlation_data = pd.concat([correlation_data, correlation_data2])

if ~{if defined(correlation_files3) then "True" else "False"}:
    correlation_filenames3 = ['~{sep="', '" correlation_files3}']
    for i in range(len(correlation_filenames3)):
        filename = correlation_filenames3[i]
        ancestry = ancestries[i]
        if ancestry not in {'AFR', 'EAS', 'NFE', 'SAS'}:
            continue
        df = pd.read_csv(filename, sep="\t", comment='#')
        if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
            print(f'{filename} has no sites in it.')
            continue
        df['ancestry'] = ancestry
        df['configuration'] = '~{configuration_label3}'
        correlation_dfs.append(df)

    correlation_data3 = pd.concat(correlation_dfs)
    correlation_data = pd.concat([correlation_data, correlation_data3])

correlation_data.to_csv('correlation_data.tsv', sep='\t')
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

ggsave(filename="correlation_plot.svg", width=7, height=5)
EOF
        Rscript script.R
    >>>

    output {
        File correlation_plot = "correlation_plot.svg"
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
		String eval_sample_id
		String truth_sample_id
		String af_expression
		String output_basename
		File af_resource
		File? sites
		String? intervals
		String? dosage_field
		Int? n_bins
		Float? right_edge_first_bin
		Float? min_af_for_accuracy_metrics
		Int mem = 16
	}
    
    Array[String] sample_map_strings = [eval_sample_id+":"+truth_sample_id]

	parameter_meta {
		evalVcf : {
			localization_optional : true
		}
		truthVcf : {
			localization_optional : true
		}
		af_resource : {
			localization_optional : true
		}
	}

	command <<<
		set -xeuo pipefail

		echo "~{truth_sample_id}:~{af_expression}" > af_expressions.list

		gatk --java-options "-Xmx~{mem - 2}G" EvaluateGenotypingPerformance -eval ~{evalVcf} -truth ~{truthVcf} --af-annotations af_expressions.list --resource ~{af_resource} \
		~{"--ids " + sites} ~{"-L " + intervals} --sample-map ~{sep=" --sample-map " sample_map_strings} ~{"--dosage-field " + dosage_field} -O ~{output_basename}.correlations.tsv \
		-OA ~{output_basename}.accuracy.tsv ~{"-nbins " + n_bins} ~{"-first-bin-right-edge " + right_edge_first_bin} ~{"--min-af-for-accuracy-metrics " + min_af_for_accuracy_metrics} --allow-differing-ploidies
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:5910defbc137d43145e6cf1f9c3539ae418b6e465250d55465e19e77773445c4"
		disks: "local-disk 100 HDD"
		memory: mem + " GB"
	}

	output {
		File correlations = "~{output_basename}.correlations.tsv"
		File accuracy = "~{output_basename}.accuracy.tsv"
	}
}