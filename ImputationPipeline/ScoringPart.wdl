version 1.0

workflow ScoringImputedDataset {
	input { 
	File weights # disease weights file. Becauase we use variant IDs with sorted alleles, there is a task at the bottom of this workflow
	#  that will allow you to sort the variants in this weights file (`SortWeights`)

	File imputed_array_vcf  # imputed VCF for scoring (and optionally PCA projection): make sure the variant IDs exactly match those in the weights file
	Int scoring_mem = 16
	Int population_scoring_mem = scoring_mem * 4
	Int vcf_to_plink_mem = 8

	String? population_basename # for naming the output of population scoring
	String basename # for naming the output of array scoring and the array projection files

	## these next 3 files are after performing PCA on your population dataset ( again, make sure all the variant IDs are the same )
	File? population_loadings
	File? population_meansd
	File? population_pcs
	File? pruning_sites_for_pca # and the sites used for PCA
	File? population_vcf # population VCF, output from PerformPopulationPCA.  The variant IDs must exactly match those in the weights file

	String? columns_for_scoring # Plink expects the first 3 columns in your weights file to be variant ID, effect allele, effect weight
	# if this isn't true, then you should give it the correct column #s in that order
	# example: if you were to set columns_for_scoring = "11 12 13" would mean that the 11th column is the variant ID, the 12th column
	# is the effect allele, and the 13th column is the effect weight
	Boolean redoPCA = false
	Boolean adjustScores = true
  }

  if (adjustScores) {
	#check for required optional inputs
	if (!defined(population_loadings)) {
		call ErrorWithMessage as ErrorPopulationLoadings {
			input:
				message = "Optional input population_loadings must be included if adjusting scores"
		}
	}

	if (!defined(population_meansd)) {
		call ErrorWithMessage as ErrorMeanSD {
			input:
				message = "Optional input population_meansd must be included if adjusting scores"
		}
	}

	if (!defined(population_pcs)) {
		call ErrorWithMessage as ErrorPopulationPCs {
			input:
				message = "Optional input population_pcs must be included if adjusting scores"
		}
	}

	if (!defined(pruning_sites_for_pca)) {
		call ErrorWithMessage as ErrorPruningSites {
			input:
				message = "Optional input pruning_sites_for_pca must be included if adjusting scores"
		}
	}

	if (!defined(population_vcf)) {
		call ErrorWithMessage as ErrorPopulationVcf {
			input:
				message = "Optional input population_vcf must be included if adjusting scores"
		}
	}
  }

	if (adjustScores) {
		call ExtractIDsPlink as ExtractIDsPopulation {
			input:
				vcf = select_first([population_vcf])
		}
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

	if (adjustScores) {
		call ExtractIDsPlink {
			input:
				vcf = imputed_array_vcf
		}

		if (redoPCA) {
			call ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
				input:
					vcf = select_first([population_vcf]),
					pruning_sites = select_first([pruning_sites_for_pca]),
					subset_to_sites = ExtractIDsPlink.ids,
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
			vcf = select_first([population_vcf]),
			basename = select_first([population_basename]),
			weights = weights,
			base_mem = population_scoring_mem,
			extra_args = columns_for_scoring,
			sites = ExtractIDsPlink.ids
		}

		call ArrayVcfToPlinkDataset {
			input:
			vcf = imputed_array_vcf,
			pruning_sites = select_first([pruning_sites_for_pca]),
			basename = basename,
			mem = vcf_to_plink_mem
		}

		call CheckPopulationIdsValid {
			input:
				pop_vcf_ids = select_first([ExtractIDsPopulation.ids]),
				pop_pc_loadings = select_first([PerformPCA.pc_loadings, population_loadings]),
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
			population_pcs = select_first([PerformPCA.pcs, population_pcs]),
			population_scores = ScorePopulation.score,
			array_pcs = ProjectArray.projections,
			array_scores = ScoreImputedArray.score
		  }
		if (!CheckPopulationIdsValid.files_are_valid) {
			call ErrorWithMessage {
				input:
				message = "Population VCF IDs are not a subset of the population PCA IDs; running with these inputs would give an incorrect result."
			}
		}
	}

  output {
	File? pc_plot = AdjustScores.pca_plot
	File? adjusted_population_scores = AdjustScores.adjusted_population_scores
	File? adjusted_array_scores = AdjustScores.adjusted_array_scores
	Boolean? fit_converged = AdjustScores.fit_converged
	File? pc_projection = ProjectArray.projections
	File raw_scores = ScoreImputedArray.score
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
	Int plink_mem = ceil(base_mem * 0.75 * 1000)
	Int disk_space =  3*ceil(size(vcf, "GB")) + 20 

	command {
		/plink2 --score ~{weights} header ignore-dup-ids list-variants-zs no-mean-imputation \
		cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --set-all-var-ids @:#:\$1:\$2 --allow-extra-chr ~{extra_args} -vcf ~{vcf} dosage=DS \
		--new-id-max-allele-len 1000 missing ~{"--extract " + sites} --out ~{basename} --memory ~{plink_mem}
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

		/plink2 --vcf ~{vcf} --extract-intersect ~{pruning_sites} ~{subset_to_sites} --allow-extra-chr --set-all-var-ids @:#:\$1:\$2 \
		--new-id-max-allele-len 1000 missing --out ~{basename} --make-bed --rm-dup force-first
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

	command <<<

		cp ~{bim} ~{basename}.bim
		cp ~{bed} ~{basename}.bed
		cp ~{fam} ~{basename}.fam

		cp ~{pc_loadings} loadings.txt
		cp ~{pc_meansd} meansd.txt

		# Check if .bim file, pc loadings, and pc meansd files have the same IDs
		# 1. extract IDs, removing first column of .bim file and first rows of the pc files
		awk '{print $2}' ~{basename}.bim > bim_ids.txt
		awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
		awk '{print $1}' meansd.txt | tail -n +2 > meansd_ids.txt

		diff bim_ids.txt pcloadings_ids.txt > diff1.txt
		diff bim_ids.txt meansd_ids.txt > diff2.txt
		diff pcloadings_ids.txt meansd_ids.txt > diff3.txt

		if [[ -s diff3.txt ]]
		then
		echo "PC loadings file and PC means file do not contain the same IDs; check your input files and run again."
		exit 1
		fi

		# check if diff files are not empty
		if [[ -s diff1.txt || -s diff2.txt ]]
		then
		echo "IDs in .bim file are not the same as the IDs in the PCA files; check that you have the right files and run again."
		exit 1
		fi

		~/flashpca/flashpca --bfile ~{basename} --project --inmeansd meansd.txt \
		--outproj projections.txt --inload loadings.txt -v
	>>>

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

			# generate the linear model from the population data using the first 4 PCs
			population_model = glm(SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = "gaussian")

			population_data$residual_score2 = resid(population_model)^2

			# generate the linear model for the variance of the score using the first 4 PCs
			population_var_model <- glm(residual_score2 ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = Gamma(link = "log"))

			# use linear model to fit full likelihood model

			# linear transformation to predict variance
			f_sigma2 <- function(t, theta) {
					PC1 = t %>% pull(PC1)
					PC2 = t %>% pull(PC2)
					PC3 = t %>% pull(PC3)
					PC4 = t %>% pull(PC4)
					PC5 = t %>% pull(PC5)
					sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
			}


			# linear transformation to predict mean
			f_mu <- function(t, theta) {
				PC1 = t %>% pull(PC1)
				PC2 = t %>% pull(PC2)
				PC3 = t %>% pull(PC3)
				PC4 = t %>% pull(PC4)
				PC5 = t %>% pull(PC5)
				mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
			}


			# negative log likelihood
			nLL_mu_and_var <- function(theta) {
				theta_mu = theta[1:5]
				theta_var = theta[6:10]
				x = population_data %>% pull(SCORE1_SUM)
				sum(log(sqrt(f_sigma2(population_data, theta_var))) + (1/2)*(x-f_mu(population_data, theta_mu))^2/f_sigma2(population_data, theta_var))
			}


			# gradient of negative log likelihood function
			grr <- function(theta) {
				theta_mu = theta[1:5]
				theta_var = theta[6:10]
				d_mu_1 <- 1
				d_mu_2 <- population_data %>% pull(PC1)
				d_mu_3 <- population_data %>% pull(PC2)
				d_mu_4 <- population_data %>% pull(PC3)
				d_mu_5 <- population_data %>% pull(PC4)
				d_sig_7 <- 1 * f_sigma2(population_data, theta_var)
				d_sig_8 <- population_data %>% pull(PC1) * f_sigma2(population_data, theta_var)
				d_sig_9 <- population_data %>% pull(PC2) * f_sigma2(population_data, theta_var)
				d_sig_10 <- population_data %>% pull(PC3) * f_sigma2(population_data, theta_var)
				d_sig_11 <- population_data %>% pull(PC4) * f_sigma2(population_data, theta_var)

				x <- population_data %>% pull(SCORE1_SUM)
				mu_coeff <- -(x - f_mu(population_data, theta_mu))/f_sigma2(population_data, theta_var)
				sig_coeff <- 1/(2*f_sigma2(population_data, theta_var)) -(1/2)*(x - f_mu(population_data, theta_mu))^2/(f_sigma2(population_data, theta_var)^2)


				grad <- c(sum(mu_coeff*d_mu_1),
				sum(mu_coeff*d_mu_2),
				sum(mu_coeff*d_mu_3),
				sum(mu_coeff*d_mu_4),
				sum(mu_coeff*d_mu_5),
				sum(sig_coeff*d_sig_7),
				sum(sig_coeff*d_sig_8),
				sum(sig_coeff*d_sig_9),
				sum(sig_coeff*d_sig_10),
				sum(sig_coeff*d_sig_11)
				)
			}

			# use linear model fits as initial parameters for full likelihood fit
			fit_mu_and_var <- optim(nLL_mu_and_var, par = c(population_model$coefficients, population_var_model$coefficients), gr = grr, method = "BFGS")


			write(ifelse(fit_mu_and_var$convergence == 0, "true", "false"), "fit_converged.txt")

		# this calculates the adjusted score for the new data
			generate_adjusted_scores = function(new_data) {
			new_data_adjusted <- new_data %>% mutate(adjusted_score = (SCORE1_SUM - f_mu(new_data, fit_mu_and_var$par[1:5]))/sqrt(f_sigma2(new_data, fit_mu_and_var$par[6:10])))
			new_data_adjusted %>% mutate(percentile=pnorm(adjusted_score,0))
			}

			# calculate adjusted score on population data,  make sure it's standardized to N(0, 1)
			population_data <- generate_adjusted_scores(population_data)



			array_scores = merge(read.csv("~{array_pcs}",  sep = "\t", header = T),
				read.csv("~{array_scores}",  sep = "\t", header = T), by.x="IID", by.y="X.IID")

			adjusted_array_scores <- generate_adjusted_scores(array_scores)

			# make sure the PCs fit well between the array and the population data
			ggplot(population_data, aes(x=PC1, y=PC2, color="Population Data")) + geom_point() + geom_point() + 
				geom_point(data = array_scores, aes(x=PC1, y=PC2, color="Array Data")) + labs(x="PC1", y="PC2") + theme_bw()
			ggsave(filename = "PCA_plot.png", dpi=300, width = 6, height = 6)

			# return population scores
			write.table(population_data %>% select(-residual_score2), file = "population_data_scores.tsv", sep="\t", row.names=F, quote = F)

			# return array scores
			write.table(adjusted_array_scores, file = "array_data_scores.tsv", sep="\t", row.names=F, quote = F)
	EOF
	>>>

	output {
		File pca_plot = "PCA_plot.png"
		File adjusted_population_scores = "population_data_scores.tsv"
		File adjusted_array_scores = "array_data_scores.tsv"
		Boolean fit_converged = read_boolean("fit_converged.txt")
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

task ExtractIDsPlink {
	input {
		File vcf
		Int disk_size = 2*ceil(size(vcf, "GB")) + 100
		Int mem = 8
	}

	Int plink_mem = ceil(mem * 0.75 * 1000)

	command <<<
		/plink2 --vcf ~{vcf} --set-all-var-ids @:#:\$1:\$2 --new-id-max-allele-len 1000 missing --write-snplist allow-dups --memory ~{plink_mem}
	>>>
	output {
		File ids = "plink2.snplist"
	}
	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem + " GB"
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

task CheckPopulationIdsValid{
	input {
		File pop_vcf_ids
		File pop_pc_loadings

	}
	command <<<
		# check if population VCF file contains a subset of population PC loading ids

		# 1. extract IDs, removing first column of .bim file and first rows of the pc files
		awk '{print $1}' ~{pop_pc_loadings} | tail -n +2 > pop__pc_ids.txt

		comm -23 <(sort pop_pc_ids.txt | uniq) <(sort ~{pop_vcf_ids} | uniq) > array_specific_ids.txt
		if [[ -s array_specific_ids.txt ]]
		then
		echo false
		else
		echo true
		fi

	>>>
	output {
		Boolean files_are_valid = read_boolean(stdout())
	}
	runtime {
		docker: "ubuntu:21.10"
	}
}

 #Print given message to stderr and return an error
 task ErrorWithMessage{
	 input {
		 String message
	 }
	command <<<
		>&2 echo "Error: ~{message}"
		exit 1
	>>>

	 runtime {
		 docker: "ubuntu:21.10"
	 }
 }


