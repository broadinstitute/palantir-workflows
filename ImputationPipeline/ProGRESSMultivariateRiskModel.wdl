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

	call ScoringTasks.ScoreVcf as ScoreImputedArray {
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
    } 

    call PCATasks.ProjectArray {
        input:
            pc_loadings = pc_loadings,
            pc_meansd = pc_meansd,
            bed = ArrayVcfToPlinkDataset.bed,
            bim = ArrayVcfToPlinkDataset.bim,
            fam = ArrayVcfToPlinkDataset.fam,
            basename = basename
    }
}

