version 1.0

struct ReferencePanelContig {
	File vcf
	File vcf_index
	File bcf
	File bcf_index
	File m3vcf
	String contig
}

struct AncestryAdjustmentModelParams {
	File fitted_model_params
	File sites_used_in_scoring
}

struct SelfExclusiveSites {
	File sites # must have columns id, chrom, pos
	Int maxAllowed
}

struct WeightSet {
	File linear_weights # standard prs weights file
	File? interaction_weights # interaction weights file, must have columns id_1, id_2, chrom_1, chrom_2, pos_1, pos_2, allele_1, allele_2, weight (order not important)
	SelfExclusiveSites? interaction_self_exclusive_sites # The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the
																												# effect alleles listed in SelfExclusizeSites.sites is observed
}

struct NamedWeightSet {
	String condition_name
	WeightSet weight_set
}

struct PRSWrapperConditionResource {
	Boolean score_condition
	Float percentile_threshold
	NamedWeightSet named_weight_set
	AncestryAdjustmentModelParams? ancestry_model_params_and_sites
}