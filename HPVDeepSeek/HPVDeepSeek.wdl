version 1.0

import "HPVDeepSeekGenotyping.wdl" as HPVDeepSeekGenotyping
import "HPVDeepSeekSomaticVariantCAlling.wdl" as HPVDeepSeekSomaticVariantCalling
import "HPVDeepSeekTertiaryAnalysis.wdl" as HPVDeepSeekTertiaryAnalysis

workflow HPVDeepSeek {
    input {
        # HPVDeepSeekGenotyping inputs
        String output_basename
        File? ubam
        File? r1_fastq
        File? r2_fastq
        File human_snp_targets_bed
        File reference
        File reference_fai
        File reference_dict
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        File capture_targets_bed
        File bait_interval_list
        String bait_set_name
        String read_group_id
        String read_group_sample_name
        String read_group_library_name = "LB_TEST"
        String read_group_platform = "ILLUMINA"
        String read_group_platform_unit = "PU_TEST"
        String read_group_description = "KAPA_TE"

        # HPVDeepSeekSomaticVariantCalling inputs
        File target_intervals
        File gnomad
        File gnomad_idx
        File pon
        File pon_idx
        File variants_for_contamination
        File variants_for_contamination_idx
        File realignment_index_bundle
        String mapping_filter_python_script = "/usr/filter_alt_ref_positions.py"
        File blastdb_nhr
        File blastdb_nin
        File blastdb_nsq
        String blastn_path = "/usr/blastn_2.2.30+"
        File funcotator_data_source
        Boolean run_alignment_artifact_filter = false

        # HPVDeepSeekTertiaryAnalysis inputs
        File high_risk_snps_hpv
        File hpv16_sublineages
    }

    call HPVDeepSeekGenotyping.HPVDeepSeekGenotyping {
        input:
            output_basename = output_basename,
            ubam = ubam,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            human_snp_targets_bed = human_snp_targets_bed,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bwa_idx_amb = bwa_idx_amb,
            bwa_idx_ann = bwa_idx_ann,
            bwa_idx_bwt = bwa_idx_bwt,
            bwa_idx_pac = bwa_idx_pac,
            bwa_idx_sa = bwa_idx_sa,
            capture_targets_bed = capture_targets_bed,
            bait_interval_list = bait_interval_list,
            bait_set_name = bait_set_name,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description
    }

    call HPVDeepSeekSomaticVariantCalling.HPVDeepSeekSomaticVariantCalling {
        input:
            output_basename = output_basename,
            tumor_bam = HPVDeepSeekGenotyping.final_bam,
            tumor_bai = HPVDeepSeekGenotyping.final_bam_index,
            target_intervals = target_intervals,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            pon = pon,
            pon_idx = pon_idx,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            realignment_index_bundle = realignment_index_bundle,
            mapping_filter_python_script = mapping_filter_python_script,
            blastdb_nhr = blastdb_nhr,
            blastdb_nin = blastdb_nin,
            blastdb_nsq = blastdb_nsq,
            blastn_path = blastn_path,
            funcotator_data_source = funcotator_data_source,
            run_alignment_artifact_filter = run_alignment_artifact_filter
    }

    call HPVDeepSeekTertiaryAnalysis.HPVDeepSeekTertiaryAnalysis {
        input:
            output_basename = output_basename,
            tumor_bam = HPVDeepSeekGenotyping.final_bam,
            tumor_bai = HPVDeepSeekGenotyping.final_bam_index,
            high_risk_snps_hpv = high_risk_snps_hpv,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            hpv16_sublineages = hpv16_sublineages
    }

    output {
        # HPVDeepSeekGenotyping outputs
        File final_bam = HPVDeepSeekGenotyping.final_bam
        File final_bam_index = HPVDeepSeekGenotyping.final_bam_index
        File umi_group_data = HPVDeepSeekGenotyping.umi_group_data
        File umi_duplication_metrics = HPVDeepSeekGenotyping.umi_duplication_metrics
        File vcf = HPVDeepSeekGenotyping.vcf
        File coverage = HPVDeepSeekGenotyping.coverage
        String top_hpv_contig = HPVDeepSeekGenotyping.top_hpv_contig
        Int top_hpv_num_reads = HPVDeepSeekGenotyping.top_hpv_num_reads
        Float top_hpv_coverage = HPVDeepSeekGenotyping.top_hpv_coverage
        Boolean is_hpv_positive = HPVDeepSeekGenotyping.is_hpv_positive
        String secondary_hpv_types = HPVDeepSeekGenotyping.secondary_hpv_types
        File fastp_report_html = HPVDeepSeekGenotyping.fastp_report_html
        File fastp_report_json = HPVDeepSeekGenotyping.fastp_report_json
        File? pre_trimmed_r1_fastqc_html = HPVDeepSeekGenotyping.r1_fastqc_html
        File? pre_trimmed_r2_fastqc_html = HPVDeepSeekGenotyping.r2_fastqc_html
        File? post_trimmed_r1_fastqc_html = HPVDeepSeekGenotyping.r1_fastqc_html
        File? post_trimmed_r2_fastqc_html = HPVDeepSeekGenotyping.r2_fastqc_html
        File pre_consensus_alignment_summary_metrics = HPVDeepSeekGenotyping.alignment_summary_metrics
        File pre_consensus_flagstat = HPVDeepSeekGenotyping.flagstat
        File pre_consensus_insert_size_metrics = HPVDeepSeekGenotyping.insert_size_metrics
        File pre_consensus_insert_size_plot = HPVDeepSeekGenotyping.insert_size_plot
        File pre_consensus_ontarget_reads = HPVDeepSeekGenotyping.ontarget_reads
        File pre_consensus_hs_metrics = HPVDeepSeekGenotyping.hs_metrics
        File pre_consensus_per_base_coverage = HPVDeepSeekGenotyping.per_base_coverage
        File post_consensus_alignment_summary_metrics = HPVDeepSeekGenotyping.alignment_summary_metrics
        File post_consensus_flagstat = HPVDeepSeekGenotyping.flagstat
        File post_consensus_insert_size_metrics = HPVDeepSeekGenotyping.insert_size_metrics
        File post_consensus_insert_size_plot = HPVDeepSeekGenotyping.insert_size_plot
        File post_consensus_ontarget_reads = HPVDeepSeekGenotyping.ontarget_reads
        File post_consensus_hs_metrics = HPVDeepSeekGenotyping.hs_metrics
        File post_consensus_per_base_coverage = HPVDeepSeekGenotyping.per_base_coverage

        # HPVDeepSeekSomaticVariantCalling outputs
        File unfiltered_vcf = HPVDeepSeekSomaticVariantCalling.unfiltered_vcf
        File unfiltered_vcf_idx = HPVDeepSeekSomaticVariantCalling.unfiltered_vcf_idx
        File mutect2_stats = HPVDeepSeekSomaticVariantCalling.mutect2_stats
        File filter_mutect_calls_stats = HPVDeepSeekSomaticVariantCalling.filtering_stats
        File filtered_vcf = HPVDeepSeekSomaticVariantCalling.output_vcf
        File filtered_vcf_idx = HPVDeepSeekSomaticVariantCalling.output_vcf_idx
        File funcotated_maf = HPVDeepSeekSomaticVariantCalling.funcotated_output_file

        # HPVDeepSeekTertiaryAnalysis outputs
        File analysis_log = HPVDeepSeekTertiaryAnalysis.analysis_log
        File breakpoints = HPVDeepSeekTertiaryAnalysis.breakpoints
        File detailed_integration_summary = HPVDeepSeekTertiaryAnalysis.detailed_integration_summary
        File integration_breakpoints = HPVDeepSeekTertiaryAnalysis.integration_breakpoints
        File integration_summary = HPVDeepSeekTertiaryAnalysis.integration_summary
        File multiple_sequence_alignment = HPVDeepSeekTertiaryAnalysis.multiple_sequence_alignment
        File phylip_formatted_msa = HPVDeepSeekTertiaryAnalysis.phylip_formatted_msa
        File phylogenetic_tree_stats = HPVDeepSeekTertiaryAnalysis.phylogenetic_tree_stats
        File phylogenetic_tree = HPVDeepSeekTertiaryAnalysis.phylogenetic_tree
        File phylogenetic_tree_visualization = HPVDeepSeekTertiaryAnalysis.phylogenetic_tree_visualization
        File sublineage_call = HPVDeepSeekTertiaryAnalysis.sublineage_call
        File high_risk_snps_found = HPVDeepSeekTertiaryAnalysis.high_risk_snps_found
    }
}