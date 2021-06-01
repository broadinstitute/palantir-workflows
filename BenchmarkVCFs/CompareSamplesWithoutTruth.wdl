version 1.0
import "FindSamplesAndBenchmark.wdl" as FindSamplesAndBenchmark

workflow CompareSamplesWithoutTruth {
  input {
    Array[String] sample_names_to_compare
    File intersected_hiconf_intervals
    Array[File] input_callset
    Array[File] ground_truth_files
    Array[File] ground_truth_indexes
    Array[String] truth_labels
    Array[File] annotation_intervals

    File gatkJarForAnnotation
    String annotationName

    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict
    File ref_fasta_sdf

    File haplotype_database
    File picard_cloud_jar

    String docker

    File? interval_list_override
    File? runs_file_override
    String comparison_docker

    String? analysis_region
  }

  Array[File] ground_truth_intervals
  scatter (i in range(length(ground_truth_files))) {
    ground_truth_intervals = flatten(ground_truth_intervals, [intersected_hiconf_intervals])
  }

  call FindSamplesAndBenchmark.FindSamplesAndBenchmark as FindSamplesAndBenchmark {
    input:
      input_callset = input_callset,
      ground_truth_files = ground_truth_files,
      ground_truth_indexes = ground_truth_indexes,
      ground_truth_intervals = ground_truth_intervals,
      truth_labels = truth_labels,
      annotation_intervals = annotation_intervals,
      skip_fingerprinting = true,
      sample_names_to_compare = sample_names_to_compare
      gatkJarForAnnotation = gatkJarFor_annotation,
      annotationName = annotationName,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_sdf = ref_fasta_sdf,
      haplotype_database = haplotype_database,
      picard_cloud_jar = picard_cloud_jar,
      docker = docker,
      interval_list_override = interval_list_override,
      runs_file_override = runs_file_override,
      comparison_docker = comparison_docker,
      analysis_region = analysis_region
  }
}
