version 1.0

import "PCARE.wdl" as PCARE
import "../Utilities/WDLs/PRSQC.wdl" as PRSQC

workflow PCAREAndQC {
    input {
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        File exome_gvcf
        File exome_gvcf_index
        File prs_weights
        Array[String]? sample_names

        File fam_history
        String basename
        
        File pc_loadings
        File pc_meansd
        File pc_sites

        Float prs_beta
        Float fam_hist_beta
        Float pc1_beta
        Float pc2_beta

        Float risk_determination_threshold_low_average = 19.69
        Float risk_determination_threshold_average_high = 20.38

        Boolean? use_ref_alt_for_ids = true
        String? chromosome_encoding = "chrMT"
        Int? mem_gb_array_vcf_to_plink = 16

        File ref_dict

        # For QC:
        File acceptable_range
        File alphashape
        Float distance_threshold
    }
    
    call PCARE.PCARE {
        input:
            imputed_wgs_vcf = imputed_wgs_vcf,
            imputed_wgs_vcf_index = imputed_wgs_vcf_index,
            exome_gvcf = exome_gvcf,
            exome_gvcf_index = exome_gvcf_index,
            prs_weights = prs_weights,
            sample_names = sample_names,
            fam_history = fam_history,
            basename = basename,
            pc_loadings = pc_loadings,
            pc_meansd = pc_meansd,
            pc_sites = pc_sites,
            prs_beta = prs_beta,
            fam_hist_beta = fam_hist_beta,
            pc1_beta = pc1_beta,
            pc2_beta = pc2_beta,
            risk_determination_threshold_low_average = risk_determination_threshold_low_average,
            risk_determination_threshold_average_high = risk_determination_threshold_average_high,
            use_ref_alt_for_ids = use_ref_alt_for_ids,
            chromosome_encoding = chromosome_encoding,
            mem_gb_array_vcf_to_plink = mem_gb_array_vcf_to_plink,
            ref_dict = ref_dict
    }

    call PRSQC.PRSQC {
        input:
            output_basename = basename,
            prs_full_risk = PCARE.full_risk,
            acceptable_range = acceptable_range,
            alphashape = alphashape,
            distance_threshold = distance_threshold
    }

    output {
        File full_risk = PCARE.full_risk
        Boolean qc_passed = PRSQC.qc_passed
        Boolean pcs_within_shape = PRSQC.pcs_within_shape
        File pca_qc_plot = PRSQC.pca_qc_plot
    }
}
