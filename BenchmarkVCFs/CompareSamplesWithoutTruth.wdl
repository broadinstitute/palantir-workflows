version 1.0
import "BenchmarkVCFs.wdl" as Benchmark
import "FindSamplesAndBenchmark.wdl" as FindSamplesAndBenchmark

workflow CompareSamplesWithoutTruth {
  input {
    Array[String] sample_names_to_compare
    File intersected_hiconf_intervals
    File input_callset
    Array[File] NYGenomes_vcf
    Array[File] NYGenomes_vcf_idx
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

    File monitoring_script

    String? analysis_region
  }

  Int VCF_disk_size = ceil(size(input_callset, "GiB") / length(sample_names_to_compare)) + 10

  # Compare samples that have truth data
  call FindSamplesAndBenchmark.FindSamplesAndBenchmark as BenchmarkFullTruthVcfs {
    input:
      input_callset = [input_callset],
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
    scatter(i in range(length(NYGenomes_vcf))) {
      call FindSamplesAndBenchmark.ExtractSampleFromCallset as ExtractFromTruth {
        input:
          callset = NYGenomes_vcf[i],
          sample = sample_name,
          basename = basename(NYGenomes_vcf[i], ".vcf.gz") + sample_name
      }
    }

    call MergeVCFs {
      input:
        input_vcfs = ExtractFromTruth.output_vcf,
        input_vcfs_indexes = ExtractFromTruth.output_vcf_index,
        output_vcf_name = "NYGC_1000G_extracted." + sample_name + ".vcf.gz"
    }

    call FindSamplesAndBenchmark.ExtractSampleFromCallset as ExtractFromInput {
      input:
        callset = input_callset,
        sample = sample_name,
        basename = basename(input_callset, ".vcf.gz") + "_extracted_" + sample_name
    }
    call FindSamplesAndBenchmark.CrosscheckFingerprints as CheckFingerprintOfExtractedSample {
      input:
        input_data = [ExtractFromInput.output_vcf],
        metrics_basename = "crosscheck",
        ground_truth_files = [MergeVCFs.output_vcf],
        haplotype_database = haplotype_database,
        disk_size = VCF_disk_size,
        preemptible_tries = 3,
        docker = docker,
        monitoring_script = monitoring_script,
        picard_jar = picard_cloud_jar
    }

    call Benchmark.Benchmark as BenchmarkVCFs {
      input:
        analysisRegion = analysis_region,
        evalVcf = ExtractFromInput.output_vcf,
        evalLabel = sample_name,
        evalVcfIndex = ExtractFromInput.output_vcf_index,
        truthVcf = MergeVCFs.output_vcf,
        confidenceInterval = intersected_hiconf_intervals,
        truthLabel = sample_name,
        truthVcfIndex = MergeVCFs.output_vcf_index,
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

    call Benchmark.Benchmark as BenchmarkDBSNP {
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

  call Benchmark.CombineSummaries as CombineSummariesWithoutTruth{
    input:
      summaries = select_all(BenchmarkVCFs.summary),
      preemptible = 1
  }

  call Benchmark.CombineSummaries as CombineSummariesDBSNP{
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

task MergeVCFs {
  input {
    File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int disk_size = 100
    Int preemptible_tries = 1
    String docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698"
    String gitc_path = "/usr/gitc/"
    Boolean no_address = false
  }
  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    bash ~{monitoring_script} > monitoring.log &
    java -Xms9000m -jar ~{gitc_path}picard.jar \
    MergeVcfs \
    INPUT=~{sep=' INPUT=' input_vcfs} \
    OUTPUT=~{output_vcf_name}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
    maxRetries: 1
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
    File monitoring_log = "monitoring.log"
  }
}
