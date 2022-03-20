version 1.0

import "PCATasks.wdl" as PCATasks
import "ScoringTasks.wdl" as ScoringTasks
import "Structs.wdl"

workflow TrainAncestryAdjustmentModel {
  input {
    NamedWeightSet named_weight_set # disease weights file

    Array[File] population_pcs # population PCs, from PerformPopulationPCA
    Array[File] population_vcfs # population VCF, output from PerformPopulationPCA.  The variant IDs must exactly match those in the weights file

    String population_basename
    File? sites # set of sites to limit scoring to
  }

  scatter (population_vcf in population_vcfs) {
    call ScoringTasks.ScoreVcf {
      input:
        vcf = population_vcf,
        basename = named_weight_set.condition_name + "_" + basename(population_vcf, ".vcf.gz"),
        weights = named_weight_set.weight_set.linear_weights,
        sites = sites
    }

    if (defined(named_weight_set.weight_set.interaction_weights)) {
      call ScoringTasks.AddInteractionTermsToScore as AddInteractionTermsToScorePopulation {
        input:
          vcf = population_vcf,
          interaction_weights = select_first([named_weight_set.weight_set.interaction_weights]),
          scores = ScoreVcf.score,
          sites = sites,
          basename = named_weight_set.condition_name + "_" + basename(population_vcf, ".vcf.gz") + population_basename,
          self_exclusive_sites = named_weight_set.weight_set.interaction_self_exclusive_sites
      }

      call ScoringTasks.CombineScoringSites {
        input:
          sites_used_linear_score = ScoreVcf.sites_scored,
          sites_used_interaction_score = AddInteractionTermsToScorePopulation.sites_used_in_interaction_score,
          basename = named_weight_set.condition_name + "_" + population_basename
      }
    }

    File scores = select_first([AddInteractionTermsToScorePopulation.scores_with_interactions, ScoreVcf.score])
    File sites_used_in_scoring_inner = select_first([CombineScoringSites.combined_scoring_sites, ScoreVcf.sites_scored])
  }

  call ScoringTasks.TrainAncestryModel {
    input:
      population_pcs = population_pcs,
      population_scores = scores,
      output_basename = named_weight_set.condition_name + "_" + population_basename
  }

  output {
    File fitted_params = TrainAncestryModel.fitted_params
    Array[File] sites_used_in_scoring = sites_used_in_scoring_inner
    File adjusted_population_scores = TrainAncestryModel.adjusted_population_scores
    Boolean fit_converged = TrainAncestryModel.fit_converged
  }
}