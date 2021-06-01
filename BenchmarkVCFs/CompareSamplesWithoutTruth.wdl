version 1.0
import "BenchmarkVCFs.wdl" as BenchmarkVCFs
import "FindSamplesAndBenchmark.wdl" as FindSamplesAndBenchmark

workflow CompareSamplesWithoutTruth {
  input {
    Array[String] sample_names_to_compare
    File intersected_hiconf_intervals
    File input_callset
    File NYGenomes_vcf
    File NYGenomes_vcf_idx
    Array[String] truth_labels
    Array[File] annotation_intervals
    Array[File] ground_truth_files
    Array[File] ground_truth_indexes
    Array [File] ground_truth_intervals
    File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
    File wgs_evaluation_regions = "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list"

    File gatkJarForAnnotation
    String annotationName

    Array[File] strat_intervals
    Array[String] strat_labels
    Array[String] jexl_variant_selectors = ["vc.isSimpleIndel()  && vc.getIndelLengths().0<0", "vc.isSimpleIndel() && vc.getIndelLengths().0>0"]
    Array[String] variant_selector_labels = ["deletion","insertion"]

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

    File monitoring_script

    String? analysis_region
  }

  Int VCF_disk_size = ceil(size(input_callset, "GiB")) + 10

  # Compare samples that have truth data
  scatter (truth_sample_name in ["NA12878","NA24143","NA24149","NA24385","NA24631","NA24694","NA24695"]) {
    call FindSamplesAndBenchmark.ExtractSampleFromCallset as ExtractTruthSamplesFromInput {
      input:
        callset = input_callset,
        sample_name = truth_sample_name,
        basename = basename(input_callset, ".vcf.gz") + "_extracted_" + truth_sample_name
    }
  }

  call FindSamplesAndBenchmark.FindSamplesAndBenchmark as BenchmarkFullTruthVcfs {
    input:
      input_callset = ExtractTruthSamplesFromInput.output_vcf,
      ground_truth_files = ground_truth_files,
      ground_truth_indexes = ground_truth_indexes,
      ground_truth_intervals = ground_truth_intervals,
      truth_labels = truth_labels,
      annotation_intervals = annotation_intervals,
      gatkJarForAnnotation = gatkJarForAnnotation,
      annotationName = annotationName,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_sdf = ref_fasta_sdf,
      haplotype_database = haplotype_database,
      picard_cloud_jar = picard_cloud_jar,
      docker = docker,
      analysis_region = analysis_region,
      stratIntervals = strat_intervals,
      stratLabels = strat_labels,
      jexlVariantSelectors = jexl_variant_selectors,
      variantSelectorLabels = variant_selector_labels,
      monitoring_script = monitoring_script
  }

  # Compare all the samples without truth data
  scatter(sample_name in sample_names_to_compare) {
    call FindSamplesAndBenchmark.ExtractSampleFromCallset as ExtractFromTruth {
      input:
        callset = NYGenomes_vcf,
        sample = sample_name,
        basename = basename(NYGenomes_vcf, ".vcf.gz") + sample_name
    }
    call FindSamplesAndBenchmark.ExtractSampleFromCallset as ExtractFromInput {
      input:
        callset = input_callset,
        sample = sample_name,
        basename = basename(input_callset, ".vcf.gz") + "_extracted_" + sample_name
    }
    call FindSamplesAndBenchmark.CrosscheckFingerprints as CheckFingerprintOfExtractedSample {
      input:
        input_data = ExtractFromInput.output_vcf,
        metrics_basename = "crosscheck",
        ground_truth_files = ExtractFromTruth.output_vcf,
        haplotype_database = haplotype_database,
        disk_size = VCF_disk_size,
        preemptible_tries = 3,
        docker = docker,
        monitoring_script = monitoring_script,
        picard_jar = picard_cloud_jar
    }

    call BenchmarkVCFs.BenchmarkVCFs as BenchmarkVCFs {
      input:
        analysisRegion = analysis_region,
        evalVcf = ExtractFromInput.output_vcf,
        evalLabel = sample_name,
        evalVcfIndex = ExtractFromInput.output_vcf_index,
        truthVcf = ExtractFromTruth.output_vcf,
        confidenceInterval = intersected_hiconf_intervals,
        truthLabel = sample_name,
        truthVcfIndex = ExtractFromTruth.output_vcf_index,
        reference = ref_fasta,
        refIndex = ref_fasta_index,
        refDict = ref_fasta_dict,
        hapMap = haplotype_database,
        stratIntervals = strat_intervals,
        stratLabels = strat_labels,
        jexlVariantSelectors = jexl_variant_selectors,
        variantSelectorLabels = variant_selector_labels,
        referenceVersion = "1",
        doIndelLengthStratification=false,
        gatkTag="4.0.11.0",
        requireMatchingGenotypes=true,
        passingOnly=true,
        vcfScoreField = "INFO.TREE_SCORE",
        gatkJarForAnnotation = gatkJarForAnnotation,
        annotationName = annotationName
    }

    call BenchmarkVCFs.BenchmarkVCFs as BenchmarkDBSNP {
      input:
        analysisRegion = analysis_region,
        evalVcf = ExtractFromInput.output_vcf,
        evalLabel = sample_name,
        evalVcfIndex = ExtractFromInput.output_vcf_index,
        truthVcf = dbsnp_vcf,
        confidenceInterval = wgs_evaluation_regions,
        truthLabel = "dbsnp",
        truthVcfIndex = dbsnp_vcf_index,
        reference = ref_fasta,
        refIndex = ref_fasta_index,
        refDict = ref_fasta_dict,
        hapMap = haplotype_database,
        stratIntervals = strat_intervals,
        stratLabels = strat_labels,
        jexlVariantSelectors = jexl_variant_selectors,
        variantSelectorLabels = variant_selector_labels,
        referenceVersion = "1",
        doIndelLengthStratification=false,
        gatkTag="4.0.11.0",
        requireMatchingGenotypes=false,
        passingOnly=true,
        vcfScoreField = "INFO.TREE_SCORE",
        gatkJarForAnnotation = gatkJarForAnnotation,
        annotationName = annotationName
    }
  }

  call BenchmarkVCFs.CombineSummaries as CombineSummariesWithoutTruth{
    input:
      summaries = select_all(BenchmarkVCFs.summary),
      preemptible = 1
  }

  call BenchmarkVCFs.CombineSummaries as CombineSummariesDBSNP{
    input:
      summaries = select_all(BenchmarkDBSNP.summary),
      preemptible = 1
  }

  output {
    File without_truth_summary = CombineSummariesWithoutTruth.summaryOut
    File dbsnp_summary = CombineSummariesDBSNP.summaryOut
    File with_truth_summary = BenchmarkFullTruthVcfs.benchmark_vcf_summary
  }
}
