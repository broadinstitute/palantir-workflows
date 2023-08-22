version 1.0
import "Structs.wdl"
import "ScoringTasks.wdl" as ScoringTasks
import "PCATasks.wdl" as PCATasks

workflow ProGRESSMultivariateRiskModel {
    input {
        File vcf
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
    }

    call ScoringTasks.DetermineChromosomeEncoding {
		input:
			weights = prs_weights
	}

	call ScoringTasks.ScoreVcf {
		input:
			vcf = vcf,
			basename = basename,
			weights = prs_weights,
			chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids
	}

    call PCATasks.ArrayVcfToPlinkDataset {
        input:
			vcf = vcf,
			pruning_sites = pc_sites,
			basename = basename,
            use_ref_alt_for_ids = use_ref_alt_for_ids,
            chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding,
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
            prs = ScoreVcf.score,
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
        
        full_risk.to_csv("~{basename}_full_risk.tsv", sep="\tsv")

        EOF

    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
    }

    output {
        File full_risk = "~{basename}_full_risk.tsv"
    }
}

