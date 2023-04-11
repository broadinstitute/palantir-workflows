version 1.0

import "CollectBGEImputationMetricsCohort.wdl" as Collect

workflow CollectBGEImputationMetricsCohortScattered {
    input {
        Array[String] contigs_hg19
        Array[String] contigs_hg38

        Array[String] ancestries
        Boolean collect_af_0_1_single_number = false

        Array[String] sample_ids1
        File eval_vcf1
        Array[String] truth_sample_ids1
        File truth_vcf1
        String configuration_label1
        File annotation_vcf1
        Map[String, String] ancestry_to_af_annotation_map1

        Array[String]? sample_ids2
        File? eval_vcf2
        Array[String]? truth_sample_ids2
        File? truth_vcf2
        String? configuration_label2
        File? annotation_vcf2
        Map[String, String] ancestry_to_af_annotation_map2

        Array[String]? sample_ids3
        File? eval_vcf3
        Array[String]? truth_sample_ids3
        File? truth_vcf3
        String? configuration_label3
        File? annotation_vcf3
        Map[String, String] ancestry_to_af_annotation_map3

        Array[String]? sample_ids4
        File? eval_vcf4
        Array[String]? truth_sample_ids4
        File? truth_vcf4
        String? configuration_label4
        File? annotation_vcf4
        Map[String, String] ancestry_to_af_annotation_map4

        Array[String]? sample_ids5
        File? eval_vcf5
        Array[String]? truth_sample_ids5
        File? truth_vcf5
        String? configuration_label5
        File? annotation_vcf5
        Map[String, String] ancestry_to_af_annotation_map5

        Array[String]? sample_ids6
        File? eval_vcf6
        Array[String]? truth_sample_ids6
        File? truth_vcf6
        String? configuration_label6
        File? annotation_vcf6
        Map[String, String] ancestry_to_af_annotation_map6

        Array[String]? sample_ids7
        File? eval_vcf7
        Array[String]? truth_sample_ids7
        File? truth_vcf7
        String? configuration_label7
        File? annotation_vcf7
        Map[String, String] ancestry_to_af_annotation_map7

        String output_basename = "cohort"
        Int plot_width = 14
        Int plot_height = 5

        Int preemptible = 1
    }

    scatter (contig in zip(contigs_hg19, contigs_hg38)) {
        call Collect.CollectBGEImputationMetricsCohort as CollectMetrics {
            input:
                ancestries = ancestries,
                collect_af_0_1_single_number = collect_af_0_1_single_number,
                output_basename = output_basename,
                plot_width = plot_width,
                plot_height = plot_height,
                preemptible = preemptible,

                sample_ids1 = sample_ids1,
                eval_vcf1 = eval_vcf1,
                truth_sample_ids1 = truth_sample_ids1,
                truth_vcf1 = truth_vcf1,
                configuration_label1 = configuration_label1,
                annotation_vcf1 = annotation_vcf1,
                ancestry_to_af_annotation_map1 = ancestry_to_af_annotation_map1,
                intervals1 = contig.right,

                sample_ids2 = sample_ids2,
                eval_vcf2 = eval_vcf2,
                truth_sample_ids2 = truth_sample_ids2,
                truth_vcf2 = truth_vcf2,
                configuration_label2 = configuration_label2,
                annotation_vcf2 = annotation_vcf2,
                ancestry_to_af_annotation_map2 = ancestry_to_af_annotation_map2,
                intervals2 = contig.right,

                sample_ids3 = sample_ids3,
                eval_vcf3 = eval_vcf3,
                truth_sample_ids3 = truth_sample_ids3,
                truth_vcf3 = truth_vcf3,
                configuration_label3 = configuration_label3,
                annotation_vcf3 = annotation_vcf3,
                ancestry_to_af_annotation_map3 = ancestry_to_af_annotation_map3,
                intervals3 = contig.right,

                sample_ids4 = sample_ids4,
                eval_vcf4 = eval_vcf4,
                truth_sample_ids4 = truth_sample_ids4,
                truth_vcf4 = truth_vcf4,
                configuration_label4 = configuration_label4,
                annotation_vcf4 = annotation_vcf4,
                ancestry_to_af_annotation_map4 = ancestry_to_af_annotation_map4,
                intervals4 = contig.left,

                sample_ids5 = sample_ids5,
                eval_vcf5 = eval_vcf5,
                truth_sample_ids5 = truth_sample_ids5,
                truth_vcf5 = truth_vcf5,
                configuration_label5 = configuration_label5,
                annotation_vcf5 = annotation_vcf5,
                ancestry_to_af_annotation_map5 = ancestry_to_af_annotation_map5,
                intervals5 = contig.right,

                sample_ids6 = sample_ids6,
                eval_vcf6 = eval_vcf6,
                truth_sample_ids6 = truth_sample_ids6,
                truth_vcf6 = truth_vcf6,
                configuration_label6 = configuration_label6,
                annotation_vcf6 = annotation_vcf6,
                ancestry_to_af_annotation_map6 = ancestry_to_af_annotation_map6,
                intervals6 = contig.right,

                sample_ids7 = sample_ids7,
                eval_vcf7 = eval_vcf7,
                truth_sample_ids7 = truth_sample_ids7,
                truth_vcf7 = truth_vcf7,
                configuration_label7 = configuration_label7,
                annotation_vcf7 = annotation_vcf7,
                ancestry_to_af_annotation_map7 = ancestry_to_af_annotation_map7,
                intervals7 = contig.right,
        }
    }
    output {
        Array[File] correlation_data = CollectMetrics.correlation_data
    }
}
