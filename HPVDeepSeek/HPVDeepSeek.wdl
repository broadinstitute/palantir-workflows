version 1.0

import "HPVDeepSeekGenotyping.wdl" as HPVDeepSeekGenotyping
import "HPVDeepSeekSomaticVariantCalling.wdl" as HPVDeepSeekSomaticVariantCalling
import "HPVDeepSeekTertiaryAnalysis.wdl" as HPVDeepSeekTertiaryAnalysis
import "HPVDeepSeekNormalization.wdl" as HPVDeepSeekNormalization

workflow HPVDeepSeek {
    input {
        # HPVDeepSeekGenotyping inputs
        String output_basename
        File r1_fastq
        File r2_fastq
        File human_snp_targets_bed
        File reference
        File reference_fai
        File reference_dict
        File bwa_idx_amb
        File bwa_idx_ann
        File bwa_idx_bwt
        File bwa_idx_pac
        File bwa_idx_sa
        File hpv_bait_interval_list
        File hpv_target_interval_list
        File hg38_bait_interval_list
        File hg38_target_interval_list
        File low_risk_hpv_genotypes
        String bait_set_name
        String read_group_id
        String read_group_sample_name
        String read_group_library_name = "LB_DEFAULT"
        String read_group_platform = "ILLUMINA"
        String read_group_platform_unit = "PU_DEFAULT"
        String read_group_description = "DS_DEFAULT"
        String read_structure = "3M2S+T"

        # HPVDeepSeekSomaticVariantCalling inputs
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

        # HPVDeepSeekNormalization inputs
        File fp_intervals
        Float ul_plasma
        Float ng_cfdna
    }

    call HPVDeepSeekGenotyping.HPVDeepSeekGenotyping { # !NameCollision
        input:
            output_basename = output_basename,
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
            hpv_bait_interval_list = hpv_bait_interval_list,
            hpv_target_interval_list = hpv_target_interval_list,
            hg38_bait_interval_list = hg38_bait_interval_list,
            hg38_target_interval_list = hg38_target_interval_list,
            low_risk_hpv_genotypes = low_risk_hpv_genotypes,
            bait_set_name = bait_set_name,
            read_group_id = read_group_id,
            read_group_sample_name = read_group_sample_name,
            read_group_library_name = read_group_library_name,
            read_group_platform = read_group_platform,
            read_group_platform_unit = read_group_platform_unit,
            read_group_description = read_group_description,
            read_structure = read_structure
    }

    call HPVDeepSeekSomaticVariantCalling.HPVDeepSeekSomaticVariantCalling { # !NameCollision
        input:
            output_basename = output_basename,
            tumor_bam = HPVDeepSeekGenotyping.duplex_bam,
            tumor_bai = HPVDeepSeekGenotyping.duplex_bam_index,
            mutect_target_intervals = hg38_target_interval_list,
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

    call HPVDeepSeekTertiaryAnalysis.HPVDeepSeekTertiaryAnalysis { # !NameCollision
        input:
            output_basename = output_basename,
            tumor_bam = HPVDeepSeekGenotyping.simplex_bam,
            tumor_bai = HPVDeepSeekGenotyping.simplex_bam_index,
            high_risk_snps_hpv = high_risk_snps_hpv,
            reference = reference,
            hpv16_sublineages = hpv16_sublineages
    }

    call HPVDeepSeekNormalization.HPVDeepSeekNormalization { # !NameCollision
        input:
            sample_id = output_basename,
            simplex_bam = HPVDeepSeekGenotyping.simplex_bam,
            simplex_bam_index = HPVDeepSeekGenotyping.simplex_bam_index,
            hpv_status = HPVDeepSeekGenotyping.hpv_status,
            fp_intervals = fp_intervals,
            ul_plasma = ul_plasma,
            ng_cfdna = ng_cfdna
    }

    output {
        # HPVDeepSeekGenotyping outputs
        File raw_bam = HPVDeepSeekGenotyping.raw_bam
        File raw_bam_index = HPVDeepSeekGenotyping.raw_bam_index
        File simplex_bam = HPVDeepSeekGenotyping.simplex_bam
        File simplex_bam_index = HPVDeepSeekGenotyping.simplex_bam_index
        File duplex_bam = HPVDeepSeekGenotyping.duplex_bam
        File duplex_bam_index = HPVDeepSeekGenotyping.duplex_bam_index
        File simplex_umi_grouped_bam = HPVDeepSeekGenotyping.simplex_umi_grouped_bam
        File simplex_umi_group_data = HPVDeepSeekGenotyping.simplex_umi_group_data
        File duplex_umi_grouped_bam = HPVDeepSeekGenotyping.duplex_umi_grouped_bam
        File duplex_umi_group_data = HPVDeepSeekGenotyping.duplex_umi_group_data
        File simplex_umi_duplication_metrics = HPVDeepSeekGenotyping.simplex_umi_duplication_metrics
        File duplex_umi_duplication_metrics = HPVDeepSeekGenotyping.duplex_umi_duplication_metrics
        File vcf = HPVDeepSeekGenotyping.vcf
        File coverage = HPVDeepSeekGenotyping.coverage
        File hpv_status = HPVDeepSeekGenotyping.hpv_status
        File fastp_report_html = HPVDeepSeekGenotyping.fastp_report_html
        File fastp_report_json = HPVDeepSeekGenotyping.fastp_report_json
        File pre_trimmed_r1_fastqc_html = HPVDeepSeekGenotyping.pre_trimmed_r1_fastqc_html
        File pre_trimmed_r2_fastqc_html = HPVDeepSeekGenotyping.pre_trimmed_r2_fastqc_html
        File post_trimmed_r1_fastqc_html = HPVDeepSeekGenotyping.post_trimmed_r1_fastqc_html
        File post_trimmed_r2_fastqc_html = HPVDeepSeekGenotyping.post_trimmed_r2_fastqc_html
        File pre_consensus_alignment_summary_metrics = HPVDeepSeekGenotyping.pre_consensus_alignment_summary_metrics
        File pre_consensus_insert_size_metrics = HPVDeepSeekGenotyping.pre_consensus_insert_size_metrics
        File pre_consensus_insert_size_plot = HPVDeepSeekGenotyping.pre_consensus_insert_size_plot
        File post_consensus_alignment_summary_metrics = HPVDeepSeekGenotyping.post_consensus_alignment_summary_metrics
        File post_consensus_insert_size_metrics = HPVDeepSeekGenotyping.post_consensus_insert_size_metrics
        File post_consensus_insert_size_plot = HPVDeepSeekGenotyping.post_consensus_insert_size_plot
        File raw_hpv_hs_metrics = HPVDeepSeekGenotyping.raw_hpv_hs_metrics
        File raw_hpv_per_target_coverage = HPVDeepSeekGenotyping.raw_hpv_per_target_coverage
        File raw_hg38_hs_metrics = HPVDeepSeekGenotyping.raw_hg38_hs_metrics
        File raw_hg38_per_target_coverage = HPVDeepSeekGenotyping.raw_hg38_per_target_coverage
        File simplex_hpv_hs_metrics = HPVDeepSeekGenotyping.simplex_hpv_hs_metrics
        File simplex_hpv_per_target_coverage = HPVDeepSeekGenotyping.simplex_hpv_per_target_coverage
        File simplex_hg38_hs_metrics = HPVDeepSeekGenotyping.simplex_hg38_hs_metrics
        File simplex_hg38_per_target_coverage = HPVDeepSeekGenotyping.simplex_hg38_per_target_coverage
        File duplex_hpv_hs_metrics = HPVDeepSeekGenotyping.duplex_hpv_hs_metrics
        File duplex_hpv_per_target_coverage = HPVDeepSeekGenotyping.duplex_hpv_per_target_coverage
        File duplex_hg38_hs_metrics = HPVDeepSeekGenotyping.duplex_hg38_hs_metrics
        File duplex_hg38_per_target_coverage = HPVDeepSeekGenotyping.duplex_hg38_per_target_coverage
        File family_sizes = HPVDeepSeekGenotyping.family_sizes
        File duplex_family_sizes = HPVDeepSeekGenotyping.duplex_family_sizes
        File duplex_yield_metrics = HPVDeepSeekGenotyping.duplex_yield_metrics
        File umi_counts = HPVDeepSeekGenotyping.umi_counts
        File duplex_qc = HPVDeepSeekGenotyping.duplex_qc

        # HPVDeepSeekSomaticVariantCalling outputs
        File contamination_table = HPVDeepSeekSomaticVariantCalling.contamination_table
        File unfiltered_vcf = HPVDeepSeekSomaticVariantCalling.unfiltered_vcf
        File unfiltered_vcf_idx = HPVDeepSeekSomaticVariantCalling.unfiltered_vcf_idx
        File mutect2_stats = HPVDeepSeekSomaticVariantCalling.mutect2_stats
        File filter_mutect_calls_stats = HPVDeepSeekSomaticVariantCalling.filter_mutect_calls_stats
        File filtered_vcf = HPVDeepSeekSomaticVariantCalling.filtered_vcf
        File filtered_vcf_idx = HPVDeepSeekSomaticVariantCalling.filtered_vcf_idx
        File funcotated_maf = HPVDeepSeekSomaticVariantCalling.funcotated_maf

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

        # HPVDeepSeekNormalization outputs
        File normalized_hpv = HPVDeepSeekNormalization.normalized_hpv
    }
}