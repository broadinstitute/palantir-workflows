version 1.0

struct GlimpseEvalSet {
  Array[String] eval_sample_ids
  Array[String] truth_sample_ids
  File eval_vcf
  File truth_vcf
  String configuration_label
  File annotation_vcf
  Map[String, String] ancestry_to_af_annoation_map1
  String? intervals
}
workflow CollectBGEImputationMetricsCohort {
  input {
    Array[String] ancestries
    Boolean collect_af_0_1_single_number = false

    Array[GlimpseEvalSet] glimpse_eval_sets

    String output_basename = "cohort"
    Int plot_width = 14
    Int plot_height = 5

    Int preemptible = 1
  }

  scatter (glimpse_eval_set in glimpse_eval_set) {
    call PearsonCorrelationByAF as PearsonByAF {
      input:
        evalVcf = glimpse_eval_set.eval_vcf,
        af_resource = glimpse_eval_set.annotation_vcf,
        ancestries = ancestries,
        ancestry_to_af_annotation_map = glimpse_eval_set.ancestry_to_af_annotation,
        truthVcf = glimpse_eval_set.truth_vcf,
        intervals = glimpse_eval_set.intervals,
        output_basename = output_basename,
        eval_sample_ids = glimpse_eval_set.eval_sample_ids,
        truth_sample_ids = glimpse_eval_set.truth_sample_ids,
        preemptible = preemptible
    }

    if (collect_af_0_1_single_number) {
      call PearsonCorrelationByAF as PearsonByAF01 {
        input:
          evalVcf = glimpse_eval_set.eval_vcf,
          af_resource = glimpse_eval_set.annotation_vcf,
          ancestries = ancestries,
          ancestry_to_af_annotation_map = glimpse_eval_set.ancestry_to_af_annotation,
          truthVcf = glimpse_eval_set.truth_vcf,
          min_af_for_accuracy_metrics = 0.1,
          n_bins = 2,
          right_edge_first_bin = 0.1,
          intervals = glimpse_eval_set.intervals,
          output_basename = output_basename,
          eval_sample_ids = glimpse_eval_set.eval_sample_ids,
          truth_sample_ids = glimpse_eval_set.truth_sample_ids,
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
    Array[File] correlation_files

    Array[String] configuration_labels
    Array[String] sample_ids
    Array[String] ancestries

    Int plot_width
    Int plot_height

    Int mem_gb = 2
    Int preemptible = 1
  }

  command <<<

    Rscript <<'EOF'
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(purrr)

    corr <- map2(c("~{sep='","' correlation_files}"), c("~{sep='","' configuration_labels}"), read_tsv(.x, skip=6) %>% mutate(configuration=.y)) %>% reduce(bind_rows)
    ancestry_map = tibble(Sample=c("~{sep="," sample_ids}"), ancestries = c("~{sep="," ancestries}"))
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
    disks: "local-disk 100 HDD"
    memory: mem_gb + " GB"
    preemptible: preemptible
  }

  output {
    File correlations = "~{output_basename}.correlations.tsv"
    File accuracy = "~{output_basename}.accuracy.tsv"
  }
}