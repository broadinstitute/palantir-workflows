version 1.0

import "cnv_germline_cohort_workflow.wdl" as cnv_germline_cohort
import "single_sample_cnv_germline_case_filter_workflow.wdl" as cnv_germline_case_filter
import "cnv_calling_and_merge_for_fabric.wdl" as cnv_calling_and_merge_for_fabric
workflow CohortCNVCallingAndMergeForFabric {
    input {
        ##################################
      #### required basic arguments ####
      ##################################
      File intervals
      Array[String]+ normal_bams
      Array[String]+ normal_bais
      Array[File]? pon_counts
      String cohort_entity_id
      File contig_ploidy_priors
      Int num_intervals_per_scatter
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker
      Array[File] short_variant_vcfs

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length

      Int? mem_gb_for_germline_cnv_caller
      Int? cpu_for_germline_cnv_caller

      Int? mem_gb_for_postprocess_germline_cnv_calls
      Int? disk_space_gb_for_postprocess_germline_cnv_calls

      Int ref_copy_number_autosomal_contigs
      

      ##########################
      #### arguments for QC ####
      ##########################
      Int maximum_number_events_per_sample
      Int maximum_number_pass_events_per_sample


      ##################################################################
      #### optional arguments for ExtractPoNFreqAnnotateFilterAndQC ####
      ##################################################################
      Array[String]? filter_expressions #= "QUAL < 50"  # Example filter criteria
      Array[String]? filter_names
      Int? mem_gb_for_extract_pon_freq
      Int? disk_for_extract_pon_freq
      Float? overlap_thresh
    }

    call cnv_germline_cohort.CNVGermlineCohortWorkflow {
        input:
            intervals = intervals,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            cohort_entity_id = cohort_entity_id,
            contig_ploidy_priors = contig_ploidy_priors,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta = ref_fasta,
            gatk_docker = gatk_docker,
            padding = padding,
            bin_length = bin_length,
            mem_gb_for_postprocess_germline_cnv_calls = mem_gb_for_postprocess_germline_cnv_calls,
            disk_space_gb_for_postprocess_germline_cnv_calls = disk_space_gb_for_postprocess_germline_cnv_calls,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,
            mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
            pon_counts = pon_counts,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs
    }

    scatter (i in range(length(CNVGermlineCohortWorkflow.genotyped_segments_vcfs))) {
        call cnv_germline_case_filter.ExtractPoNFreqAnnotateFilterAndQC {
            input:
                vcf = CNVGermlineCohortWorkflow.genotyped_segments_vcfs[i],
                vcf_idx = CNVGermlineCohortWorkflow.genotyped_segments_vcf_indexes[i],
                panel_vcfs = select_first([CNVGermlineCohortWorkflow.pon_genotyped_segments_vcfs, CNVGermlineCohortWorkflow.genotyped_segments_vcfs]),
                intervals = CNVGermlineCohortWorkflow.preprocessed_intervals,
                overlap_thresh = overlap_thresh,
                filter_expressions = filter_expressions,
                filter_names = filter_names,
                max_events = maximum_number_events_per_sample,
                max_pass_events = maximum_number_pass_events_per_sample,
                sample_name = CNVGermlineCohortWorkflow.read_counts_entity_ids[i],
                mem_gb = mem_gb_for_extract_pon_freq,
                disk_size_gb = disk_for_extract_pon_freq
        }

        Boolean qc_passed_scatter = (ExtractPoNFreqAnnotateFilterAndQC.qc_status == "PASS")

        call cnv_calling_and_merge_for_fabric.ReformatAndMergeForFabric {
            input:
                cnv_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf,
                short_variant_vcf = short_variant_vcfs[i],
                gatk_docker = gatk_docker
        }

        call cnv_calling_and_merge_for_fabric.GCNVVisualzation {
            input:
            filtered_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf,
            case_copy_ratios = CNVGermlineCohortWorkflow.denoised_copy_ratios[i],
            case_read_counts = CNVGermlineCohortWorkflow.read_counts[i],
            panel_copy_ratios = select_first([CNVGermlineCohortWorkflow.pon_denoised_copy_ratios, CNVGermlineCohortWorkflow.denoised_copy_ratios]),
            panel_read_counts = pon_counts,
            interval_lists = CNVGermlineCohortWorkflow.sharded_interval_lists
        }
    }

    output {
        Array[File] filtered_cnv_genotyped_segments_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf
        Array[File] filtered_cnv_genotyped_segments_vcf_index = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_index
        Array[File] filtered_cnv_genotyped_segments_vcf_md5sum = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_md5sum

        Array[File] merged_vcf = ReformatAndMergeForFabric.merged_vcf
        Array[File] merged_vcf_index = ReformatAndMergeForFabric.merged_vcf_index
        Array[File] merged_vcf_md5sum = ReformatAndMergeForFabric.merged_vcf_md5sum

        Array[Boolean] qc_passed = qc_passed_scatter
        Array[File] cnv_metrics = ExtractPoNFreqAnnotateFilterAndQC.cnv_metrics
        Array[File] cnv_event_report = GCNVVisualzation.cnv_event_report

    }

}
