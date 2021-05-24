version 1.0

workflow ScoringImputedDataset {
	input { 
    File weights # disease weights file. Becauase we use variant IDs with sorted alleles, there is a task at the bottom of this workflow
    #  that will allow you to sort the variants in this weights file (`SortWeights`)

    File imputed_array_vcf  # imputed VCF for scoring (and optionally PCA projection): make sure the variant IDs exactly match those in the weights file
    Int scoring_mem = 16
    Int population_scoring_mem = scoring_mem * 4
    Int vcf_to_plink_mem = 8
    
    String population_basename # for naming the output of population scoring
    String basename # for naming the output of array scoring and the array projection files

    ## these next 3 files are after performing PCA on your population dataset ( again, make sure all the variant IDs are the same )
    File population_loadings
    File population_meansd
    File population_pcs
    File pruning_sites_for_pca # and the sites used for PCA
    File population_vcf # population VCF, output from PerformPopulationPCA.  The variant IDs must exactly match those in the weights file

    File adjust_scores_rscript = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/ScoringAdjustment.R" 
    String? columns_for_scoring # Plink expects the first 3 columns in your weights file to be variant ID, effect allele, effect weight
    # if this isn't true, then you should give it the correct column #s in that order
    # example: if you were to set columns_for_scoring = "11 12 13" would mean that the 11th column is the variant ID, the 12th column 
    # is the effect allele, and the 13th column is the effect weight
    Boolean redoPCA = false
  }

  call ExtractIDs as ExtractIDsPopulation {
  	input:
  		vcf = population_vcf,
  		output_basename = population_basename
  }
  
  call ScoreVcf as ScoreImputedArray {
  	input:
  	vcf = imputed_array_vcf,
  	basename = basename,
  	weights = weights,
  	base_mem = scoring_mem,
  	extra_args = columns_for_scoring,
  	sites = ExtractIDsPopulation.ids
  }

  call ExtractIDs {
  	input:
  		vcf = imputed_array_vcf,
  		output_basename = basename
  }

  if (redoPCA) {
  	call ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
  		input:
  			vcf = population_vcf,
  			pruning_sites = pruning_sites_for_pca,
  			subset_to_sites = ExtractIDs.ids,
  			basename = "population"
  	}

  	call PerformPCA {
        input:
          bim = PopulationArrayVcfToPlinkDataset.bim,
          bed = PopulationArrayVcfToPlinkDataset.bed,
          fam = PopulationArrayVcfToPlinkDataset.fam,
          basename = basename
      }
  }

  call ScoreVcf as ScorePopulation {
  	input:
  	vcf = population_vcf,
  	basename = population_basename,
  	weights = weights,
  	base_mem = population_scoring_mem,
  	extra_args = columns_for_scoring,
  	sites = ExtractIDs.ids
  }

  call ArrayVcfToPlinkDataset {
  	input:
  	vcf = imputed_array_vcf,
  	pruning_sites = pruning_sites_for_pca,
  	basename = basename,
  	mem = vcf_to_plink_mem
  }

  call ProjectArray {
  	input:
  	pc_loadings = select_first([PerformPCA.pc_loadings, population_loadings]),
  	pc_meansd = select_first([PerformPCA.mean_sd, population_meansd]),
  	bed = ArrayVcfToPlinkDataset.bed,
  	bim = ArrayVcfToPlinkDataset.bim,
  	fam = ArrayVcfToPlinkDataset.fam,
  	basename = basename
  }

  call AdjustScores {
  	input:
  	adjusting_Rscript = adjust_scores_rscript,
  	population_pcs = select_first([PerformPCA.pcs, population_pcs]),
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
		File? sites
	}

	Int runtime_mem = base_mem + 2
	Int plink_mem = base_mem * 1000
	Int disk_space =  3*ceil(size(vcf, "GB")) + 20 

	command {
		/plink2 --score ~{weights} header ignore-dup-ids list-variants-zs no-mean-imputation \
		cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --allow-extra-chr ~{extra_args} -vcf ~{vcf} dosage=DS \
		~{"--extract " + sites} --out ~{basename} --memory ~{plink_mem}
	}

	output {
		File score = "~{basename}.sscore"
		File log = "~{basename}.log"
		File sites_scored = "~{basename}.sscore.vars.zst"
	}

	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk " + disk_space + " HDD"
		memory: runtime_mem + " GB"
	}
}

# This just turns the array vcf into bim/bed/fam format and extracts only the sites
# Used for the original PCA steps
task ArrayVcfToPlinkDataset {
	input {
		File vcf
		File pruning_sites
		File? subset_to_sites
		String basename
		Int mem = 8
	}

	Int disk_space =  3*ceil(size(vcf, "GB")) + 20 

	command {

		/plink2 --vcf ~{vcf} --extract-intersect ~{pruning_sites} ~{subset_to_sites} --allow-extra-chr \
		--out ~{basename} --make-bed --rm-dup force-first
	}

	output {
		File bed = "~{basename}.bed"
		File bim = "~{basename}.bim"
		File fam = "~{basename}.fam"
	}

	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk " + disk_space + " HDD"
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
		docker: "quay.io/ckachuli/flashpca@sha256:85e9ee91bc552e46a0d69cc851b893419c8de6588c696458fc770eee526e381d" # a special version of flashpca which allows to project a single sample without erroring out at an unnecessary check
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

	command <<<
		Rscript - <<- "EOF"
			library(ggplot2)
			library(dplyr)

			population_pcs = read.csv("~{population_pcs}", sep="\t", header = T)
			population_scores = read.csv("~{population_scores}", sep="\t", header = T)

			population_data = merge(population_pcs, population_scores, by.x="IID", by.y="X.IID")

			# generate the linear model from the population data using the first 5 PCs
			population_model = glm(SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4 + PC5, data = population_data, family = "gaussian")

			population_data$residual_score = resid(population_model)
			population_resid_mean = mean(population_data$residual_score)
			population_resid_sd = sd(population_data$residual_score)

			# calculate adjusted score on population data,  make sure it's standardized to N(0, 1)
			population_data$adjusted_score = (population_data$residual_score - population_resid_mean)/population_resid_sd

			# this calculates the adjusted score for the new data
			generate_adjusted_scores = function(new_data) {
			subset_data_for_model = new_data %>% transmute(raw_score = SCORE1_SUM, PC1=PC1, PC2=PC2, PC3=PC3, PC4=PC4, PC5=PC5)
			new_data$residual_score = subset_data_for_model$raw_score - predict(population_model, subset_data_for_model) # calculate the residual score from the model
			new_data$adjusted_score = (new_data$residual_score - population_resid_mean)/population_resid_sd # again adjust compared to population
			new_data %>% rowwise() %>% mutate(percentile=pnorm(adjusted_score, round(mean(population_data$adjusted_score),5),
			round(sd(population_data$adjusted_score), 5) == 1))
			}

			array_scores = merge(read.csv("~{array_pcs}",  sep = "\t", header = T),
				read.csv("~{array_scores}",  sep = "\t", header = T), by.x="IID", by.y="X.IID")
	>>>

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

 task ExtractIDs {
     input {
         File vcf
         String output_basename
         Int disk_size = 2*ceil(size(vcf, "GB")) + 100
     }

     command <<<
         bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.original_array.ids
     >>>
     output {
         File ids = "~{output_basename}.original_array.ids"
     }
     runtime {
         docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
         disks: "local-disk " + disk_size + " HDD"
         memory: "4 GB"
     }
 }

 task PerformPCA {
   input {
     File bim
     File bed
     File fam
     String basename
     Int mem = 8
   }

   # again, based on Wallace commands
   command {
     cp ~{bim} ~{basename}.bim
     cp ~{bed} ~{basename}.bed
     cp ~{fam} ~{basename}.fam

     ~/flashpca/flashpca --bfile ~{basename} -n 16 -d 20 --outpc ${basename}.pc \
     --outpve ${basename}.pc.variance --outload ${basename}.pc.loadings \
     --outmeansd ${basename}.pc.meansd
   }

   output {
     File pcs = "${basename}.pc"
     File pc_variance = "${basename}.pc.variance"
     File pc_loadings = "${basename}.pc.loadings"
     File mean_sd = "${basename}.pc.meansd"
     File eigenvectors = "eigenvectors.txt"
     File eigenvalues = "eigenvalues.txt"
   }

   runtime {
     docker: "skwalker/flashpca:v1"
     disks: "local-disk 400 HDD"
     memory: mem + " GB"
   }
 }


