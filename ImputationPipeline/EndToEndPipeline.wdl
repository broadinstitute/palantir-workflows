version 1.0

import "Imputation.wdl" as imputation_pipeline
import "ScoringPart.wdl" as scoring_pipeline
import "PerformPopulationPCA.wdl" as population_pipeline # Like Thousand Genomes
import "Structs.wdl" as Structs

workflow EndToEndPipeline {
	input {
	  File multi_sample_vcf ## the dataset you want to impute, this is used both in imputation as well 
	  # as for doing the scoring adjustment steps, and may be used to select sites before LD pruning
	  # if generating new population PCs
	  File multi_sample_vcf_index ## index for the dataset you want to impute,used for imputation
	  
	  Boolean perform_extra_qc_steps # these are optional additional extra QC steps that can be performed
	  # on the `multi_sample_vcf` before imputing -- these should only be run for large and disparate sample 
	  # sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
	 

	  ## Do you want to generate new PCs for your population, either if: 
	  ## you have a new population (don't want to run Thousand Genomes)
	  ## or you have a new array and want to limit the projection to sites called in that array
	  Boolean generate_population_pcs = false


	  # the population vcf to use for the scoring adjustment, and if `generate_population_pcs` is set to
	  # true, for generating the population pcs on
	  # this defaults to Thousand Genomes
	  File population_vcf =  "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/thousand_genomes_sorted_variant_ids.vcf.gz"
	  File population_vcf_index =  "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/thousand_genomes_sorted_variant_ids.vcf.gz.tbi"

	  ## if generate_population_pcs is not true: then you need to include the PC files in order to project the
	  ## array data onto that population -- here we default to running on the Thousand Genomes with the sites
	  ## Wallace generated after LD pruning

	  File? population_loadings = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/WallacesPCASites/sorted_thousand_genomes_wallace_sites.pc.loadings"
	  File? population_meansd = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/WallacesPCASites/sorted_thousand_genomes_wallace_sites.pc.meansd"
	  File? population_pcs = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/WallacesPCASites/sorted_thousand_genomes_wallace_sites.pc"
	  File? pruning_sites_for_pca = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/RiskScoreAdjustmentFiles/WallacesPCASites/wallace_pruning_sites_sorted_ids.txt"

	  ## The following are inputs for scoring and performing the adjustment

	  WeightSet disease_weights  # disease weights file. Because we use variant IDs with sorted alleles, there is a task at the bottom of this workflow
	  String? columns_for_scoring # Plink expects the first 3 columns in your weights file to be variant ID, effect allele, effect weight
	  
	  Int scoring_mem = 16 # update memory for scoring imputed array
	  Int population_scoring_mem = 16 # update memory for scoring population VCF

	  # output names (what the files will be named)
	  String output_callset_name # the name for the imputed callset name
	  String population_basename  # the basename for the output PCs if `generate_population_pcs` is true
	  Boolean adjustScores = true

	}

  call imputation_pipeline.ImputationPipeline as ImputationSteps {
  	input:
  	  multi_sample_vcf = multi_sample_vcf,
  	  multi_sample_vcf_index = multi_sample_vcf_index,
  	  perform_extra_qc_steps = perform_extra_qc_steps,
  	  output_callset_name = output_callset_name,
  }

  if (generate_population_pcs) {
  	call population_pipeline.PerformPopulationPCA as PopulationPCASteps {
  	  input:  
  	    population_vcf = population_vcf,
  	    population_vcf_index = population_vcf_index,
  	    basename = population_basename,
  	    imputed_array_vcfs =  [ImputationSteps.imputed_multisample_vcf],
  	    original_array_vcfs = [multi_sample_vcf]
    }
  }

  call scoring_pipeline.ScoringImputedDataset as ScoringSteps {
	  input :
		weight_set  = disease_weights,
		columns_for_scoring = columns_for_scoring,
		imputed_array_vcf = ImputationSteps.imputed_multisample_vcf,
	    scoring_mem = scoring_mem,
	    population_scoring_mem = population_scoring_mem,
	    population_vcf = select_first([PopulationPCASteps.sorted_variant_id_dataset, population_vcf]),
	    population_basename = population_basename,
	    basename = output_callset_name,
	    population_loadings = select_first([PopulationPCASteps.population_loadings, population_loadings]), # either use your newely generated PC files or the input loadings/meansd/pcs/pruning sites
	    population_meansd = select_first([PopulationPCASteps.population_meansd, population_meansd]), 
	    population_pcs = select_first([PopulationPCASteps.population_pcs, population_pcs]),
	    pruning_sites_for_pca = select_first([PopulationPCASteps.pruning_sites_for_pca, pruning_sites_for_pca]),
	    adjustScores = adjustScores
	}

	output {

		# the final imputed VCF + index
		File imputed_multisample_vcf = ImputationSteps.imputed_multisample_vcf
    File imputed_multisample_vcf_index = ImputationSteps.imputed_multisample_vcf_index

    # the following files will only be generated if `generate_population_pcs` is set to true
    # these can be used in future runs to just calculate the scores and the scoring adjustment
    File? new_population_loadings = PopulationPCASteps.population_loadings
    File? new_population_meansd = PopulationPCASteps.population_loadings
    File? new_population_pcs = PopulationPCASteps.population_loadings
    File? new_pruning_sites_for_pca = PopulationPCASteps.population_loadings
    File? new_population_dataset_with_sorted_variant_ids = PopulationPCASteps.sorted_variant_id_dataset
    File? new_population_dataset_index_with_sorted_variant_ids = PopulationPCASteps.sorted_variant_id_dataset_index
    
		# results from the scoring part
    File? pc_plot = ScoringSteps.pc_plot
  	File? adjusted_population_scores = ScoringSteps.adjusted_population_scores
  	File? adjusted_array_scores = ScoringSteps.adjusted_array_scores
  	File raw_scores = ScoringSteps.raw_scores
	}
}