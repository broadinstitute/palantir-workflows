version 1.0

struct ReferencePanelContig {
	File vcf
	File vcf_index
	File bcf
	File bcf_index
	File m3vcf
	String contig
}


struct SelfExclusiveSites {
	File sites
	Int maxAllowed
}

struct WeightSet {
	File linear_weights # standard prs weights file
	File? interaction_weights # interaction weights file
	SelfExclusiveSites? interaction_self_exclusive_sites # The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the
																												# effect alleles listed in SelfExclusizeSites.sites is observed
}

struct NamedWeightSet {
	String condition_name
	WeightSet weight_set
}