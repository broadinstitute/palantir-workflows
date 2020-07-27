version 1.0

workflow ScoringImputedDataset {
	input { 
    File weights # disease weights file. Becauase we use variant IDs with sorted alleles, there is a task at the bottom of this workflow
    #  that will allow you to sort the variants in this weights file

    File? original_array_vcf # original array for better PCA projection ( this is optional, with the default being it will run on the imputed array vcf instead)
    File imputed_array_vcf  # imputed VCF for scoring (and optionally PCA projection): make sure the variant IDs exactly match those in the weights file 
    Int scoring_mem = 16
    File population_vcf # population VCF (e.g., Thousand Genomes): again, make sure the variant IDs exactly match those in the weights file 
    
    String population_basename # for naming the output of population scoring
    String basename # for naming the output of array scoring and the array projection files

    ## these next 3 files are after performing PCA on your population dataset ( again, make sure all the variant IDs are the same )
    File population_loadings
    File population_meansd
    File population_pcs
    File pruning_sites_for_pca # and the sites used for PCA

    File adjust_scores_rscript = "gs://broad-dsde-methods-skwalker/ScoringAdjustment.R" 
    String? columns_for_scoring # Plink expects the first 3 columns in your weights file to be variant ID, effect allele, effect weight
    # if this isn't true, then you should give it the correct column #s in that order
    # example: if you were to set columns_for_scoring = "11 12 13" would mean that the 11th column is the variant ID, the 12th column 
    # is the effect allele, and the 13th column is the effect weight
  }

  # this adds in the correct IDs (sorted) so we can run things properly
  if (defined(original_array_vcf)) {
  	call UpdateVariantIds {
  		input:
  		  vcf = original_array_vcf, 
  		  basename = "original_array.different_ids"
  	}

  	call SortIds as SortOriginalArrayVariantIds {
  		input:
  		  vcf = UpdateVariantIds.output_vcf,
  		  basename = "original_array.updated_ids"
  	}
  }
  
  call ScoreVcf as ScoreImputedArray {
  	input:
	  	vcf = imputed_array_vcf,
	  	basename = basename,
	  	weights = weights,
	  	base_mem = scoring_mem,
	  	extra_args = columns_for_scoring
  }

  call ScoreVcf as ScorePopulation {
  	input:
  	vcf = population_vcf,
  	basename = population_basename,
  	weights = weights,
  	base_mem = scoring_mem * 4,
  	extra_args = columns_for_scoring 
  }

  call ArrayVcfToPlinkDataset {
  	input:
  	vcf = select_first([SortOriginalArrayVariantIds.output_vcf, imputed_array_vcf]),
  	pruning_sites = pruning_sites_for_pca,
  	basename = basename
  }

  call ProjectArray {
  	input:
  	pc_loadings = population_loadings,
  	pc_meansd = population_meansd,
  	bed = ArrayVcfToPlinkDataset.bed,
  	bim = ArrayVcfToPlinkDataset.bim,
  	fam = ArrayVcfToPlinkDataset.fam,
  	basename = basename
  }

  call AdjustScores {
  	input:
  	adjusting_Rscript = adjust_scores_rscript,
  	population_pcs = population_pcs,
  	population_scores = ScorePopulation.score,
  	array_pcs = ProjectArray.projections,
  	array_scores = ScoreImputedArray.score
  }

  output {
  	File pc_plot = AdjustScores.pca_plot
  	File adjusted_population_scores = AdjustScores.adjusted_population_scores
  	File adjusted_array_scores = AdjustScores.adjusted_array_scores
  }
}

# Wallace suggested using variant IDs as chrom:pos:allele1:allele2 where allele1 is < allele2 (alphabetically sorted) so 
# that the variant IDs exactly match in the disease weights file and the imputed VCF file (assuming your weights file also
# contains variant IDs sorted this way)
# If you're using the partners VCF I imputed, those variants have already been sorted in this order, so this task isn't 
# necessary (which is why it's not actually run in the workflow)

## This is assuming that your variant IDs are in the format chr:positionr:allele1:allele2 and just alphabetizes allele1 and 
## allele 2. To get variant IDs in this format, you can run `bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ~{vcf}`
task SortIds {
	input {
		File vcf
		String basename
		Int disk_space =  3*ceil(size(vcf, "GB"))
	}

	command <<<
	zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
	>>>

	output {
		File output_vcf = "~{basename}.vcf.gz"
	}

	runtime {
		docker: "skwalker/imputation:with_vcftools" 
		disks: "local-disk " + disk_space + " HDD"
		memory: "16 GB"
	}
}

# score with plink2
task ScoreVcf {
	input {
		File vcf
		String basename
		File weights
		Int base_mem = 8
		String? extra_args
	}

	Int runtime_mem = base_mem + 2
	Int plink_mem = base_mem * 1000



	command {
		/plink2 --score ~{weights} header ignore-dup-ids list-variants-zs no-mean-imputation \
		cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums -vcf ~{vcf} dosage=DS \
		--out ~{basename} --memory ~{plink_mem} ~{extra_args}
	}

	output {
		File score = "~{basename}.sscore"
		Array[File] outputs = glob("~{basename}*")
	}

	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk 400 HDD"
		memory: runtime_mem + " GB"
	}
}

# This just turns the array vcf into bim/bed/fam format and extracts only the sites
# Used for the original PCA steps
task ArrayVcfToPlinkDataset {
	input {
		File vcf
		File pruning_sites
		String basename
		Int mem = 8
	}


	command {

		/plink2 --vcf ~{vcf} --extract ~{pruning_sites} \
		--out ~{basename} --make-bed
	}

	output {
		File bed = "~{basename}.bed"
		File bim = "~{basename}.bim"
		File fam = "~{basename}.fam"
	}

	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk 400 HDD"
		memory: mem + " GB"
	}
} 

# This projects the array dataset using the previously generated PCs, using flashPCA
task ProjectArray {
	input {
		File bim
		File bed
		File fam
		File pc_loadings
		File pc_meansd
		String basename
		Int mem = 8
	}

	command {

		cp ~{bim} ~{basename}.bim
		cp ~{bed} ~{basename}.bed
		cp ~{fam} ~{basename}.fam

		cp ~{pc_loadings} loadings.txt
		cp ~{pc_meansd} meansd.txt

		~/flashpca/flashpca --bfile ~{basename} --project --inmeansd meansd.txt \
		--outproj projections.txt --inload loadings.txt -v
	}

	output {
		File projections = "projections.txt"
	}

	runtime {
		docker: "skwalker/flashpca:v1"
		disks: "local-disk 400 HDD"
		memory: mem + " GB"
	} 
}

# This does the scoring adjustment
task AdjustScores {
	input {
		File adjusting_Rscript
		File population_pcs
		File population_scores 
		File array_pcs
		File array_scores
		Int mem = 2
	}

	command {
		Rscript ~{adjusting_Rscript} ~{population_pcs} ~{population_scores} ~{array_pcs} ~{array_scores}
	}

	output {
		File pca_plot = "PCA_plot.png"
		File adjusted_population_scores = "population_data_scores.tsv"
		File adjusted_array_scores = "array_data_scores.tsv"
	}

	runtime {
		docker: "skwalker/rscripting:with_rutils"
		disks: "local-disk 400 HDD"
		memory: mem + " GB"
	}
} 


# This task allows you to sort each variant ID in your weights file. It already assumes they are in the format chr:position:a1:a2
# and just sorts a1, a2. You will have to perform other awk magic to get it into this format otherwise.

task SortWeights {
  input {
    File weights_file
    Int disk_space = 50
    Int id_column # the column # of the variant IDs
     String basename # what you wanted the new weights file to be called
   }
    
    command <<<
  
    awk -v id_col="~{id_column}" -v OFS='\t' '{split($id_col, n, ":"); if ( n[4] < n[3])  $id_col=n[1]":"n[2]":"n[4]":"n[3]; print $0}' ~{weights_file} > ~{basename}.txt

	>>>

	output {
		File sorted_weights = "~{basename}.txt"
	}

	runtime {
		docker: "skwalker/imputation:with_vcftools" 
		disks: "local-disk " + disk_space + " HDD"
		memory: "16 GB"
	}
}

task UpdateVariantIds {

  input {
    File? vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~{vcf} -O z -o ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

