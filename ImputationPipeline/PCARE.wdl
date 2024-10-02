version 1.0

import "Structs.wdl"
import "PCATasks.wdl" as PCATasks
import "ScoreBGE/ScoreBGE.wdl" as ScoreBGE

workflow PCARE {
    input {
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        File exome_gvcf
        File exome_gvcf_index
        File prs_weights

        File fam_history
        String basename
        
        File pc_loadings
        File pc_meansd
        File pc_sites

        Float prs_beta
        Float fam_hist_beta
        Float pc1_beta
        Float pc2_beta

        Boolean use_ref_alt_for_ids = true
        String chromosome_encoding = "chrMT"
        Int mem_gb_array_vcf_to_plink = 16

        File ref_dict
    }

    call ScoreBGE.ScoreBGE {
        input:
            exome_gvcf = exome_gvcf,
            exome_gvcf_index = exome_gvcf_index,
            imputed_wgs_vcf = imputed_wgs_vcf,
            imputed_wgs_vcf_index = imputed_wgs_vcf_index,
            basename = basename,
            weights = prs_weights,
            score_haploid_as_diploid = true,
            ref_dict = ref_dict
    }


    call PCATasks.ArrayVcfToPlinkDataset {
        input:
            vcf = imputed_wgs_vcf,
            pruning_sites = pc_sites,
            basename = basename,
            use_ref_alt_for_ids = use_ref_alt_for_ids,
            chromosome_encoding = chromosome_encoding,
            mem = mem_gb_array_vcf_to_plink
    }

    call PCATasks.ProjectArray {
        input:
            pc_loadings = pc_loadings,
            pc_meansd = pc_meansd,
            bed = ArrayVcfToPlinkDataset.bed,
            bim = ArrayVcfToPlinkDataset.bim,
            fam = ArrayVcfToPlinkDataset.fam,
            basename = basename,
            divisor = "none"
    }

    call ComputeRiskValue {
        input:
            prs = ScoreBGE.score,
            pcs = ProjectArray.projections,
            family_history = fam_history,
            prs_beta = prs_beta,
            fam_hist_beta = fam_hist_beta,
            pc1_beta = pc1_beta,
            pc2_beta = pc2_beta,
            basename = basename
    }

    output {
        File full_risk = ComputeRiskValue.full_risk
    }
}

task ComputeRiskValue {
    input {
        File prs
        File pcs
        File family_history

        Float prs_beta
        Float fam_hist_beta
        Float pc1_beta
        Float pc2_beta

        String basename
    }

    command <<<
        python3 << "EOF"
        import pandas as pd
        prs = pd.read_csv("~{prs}", sep="\t")
        pcs = pd.read_csv("~{pcs}", sep="\t").set_index("IID")
        fam_hist = pd.read_csv("~{family_history}", sep="\t").set_index("sample_id")

        full_risk = prs.join(pcs, on='#IID').join(fam_hist, on='#IID')

        full_risk['combined_risk_score'] = (~{prs_beta}*full_risk.SCORE1_SUM + ~{fam_hist_beta}*full_risk.fam_hist + 
                                            ~{pc1_beta}*full_risk.PC1 + ~{pc2_beta}*full_risk.PC2)

        full_risk['risk_determination'] = full_risk['combined_risk_score'].apply(
            lambda score: 'low' if score < 19.69 else ('high' if score > 20.38 else 'average'))
        
        full_risk = full_risk.rename(columns={"#IID": "sample_id", "SCORE1_SUM": "prs_score", "PC1": "pc1", "PC2": "pc2", "fam_hist": "family_history"})

        full_risk.to_csv("~{basename}_full_risk.tsv", sep="\t", index=False)

        EOF

    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
    }

    output {
        File full_risk = "~{basename}_full_risk.tsv"
    }
}
