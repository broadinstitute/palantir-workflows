version 1.0

import "ScoringTasks.wdl" as ScoringTasks
import "Structs.wdl"

workflow PolygenicRiskScore {
	input { 
	File linear_weights

	File vcf  # VCF for scoring
	Int scoring_mem = 16
	Int vcf_to_plink_mem = 8

	String basename # for naming the output of array scoring and the array projection files

    Boolean assume_missing_hom_ref = false # this can be used when using a whole genome vcf, where any uncalled sites can be assumed to be hom-ref.
	#  In this case, you must also provide ref_fasta and ref_fai
	File? ref_fasta
    File? ref_fai
  }
	call ScoringTasks.DetermineChromosomeEncoding {
		input:
			weights = linear_weights
	}

	call ScoringTasks.ScoreVcf as ScoreImputedArray {
		input:
			vcf = vcf,
			basename = basename,
			weights = linear_weights,
			base_mem = scoring_mem,
			chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding,
			assume_missing_hom_ref = assume_missing_hom_ref,
			ref_fasta = ref_fasta,
			ref_fai = ref_fai
	}

  output {
	File raw_scores = ScoreImputedArray.score
  }
}