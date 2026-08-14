# PRS Scoring Validation

Workflows used to regression-test changes to the PRS scoring pipeline. See the [PRS README](../README.md) for the
pipeline these validate.

* [ValidateScoring](#validatescoring)
* [SubsetWeightSet](#subsetweightset)

## ValidateScoring

Defined in [ValidateScoring.wdl](ValidateScoring.wdl); registered on Dockstore as **ValidateScoring** (see
[../../.dockstore.yml](../../.dockstore.yml)).

### Summary

Runs the scoring pipeline several different ways over the same inputs and compares the resulting adjusted scores. It
imports `ScoringPart.wdl` twice: once from the local checkout (the branch under test) and once directly from
`https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/PRS/ScoringPart.wdl` (the `main` branch), so
one run compares a branch against `main`.

For each `NamedWeightSet` in `named_weight_sets` the workflow:

- extracts the variant IDs present in the imputed array VCF and subsets the weight set to them
  ([SubsetWeightSet](#subsetweightset));
- removes any site with a no-call genotype from the WGS VCF (`QCSites`, `vcftools --max-missing-count 0`) so that
  no-calls cannot skew the WGS score in a way the ancestry adjustment would not account for, and subsets the weight
  set again to the surviving sites;
- re-runs the population PCA over only the sites present in the imputed VCF (`PCATasks.ArrayVcfToPlinkDataset` +
  `PCATasks.PerformPCA`), since the supplied population PCA may not be appropriate for this site set;
- scores the imputed array on the branch (full weight set, and again with the WGS-restricted weight set), scores the
  imputed array on `main`, and scores the WGS VCF with the WGS-restricted weight set;
- separately trains the ancestry-adjustment model with `TrainAncestryAdjustmentModel` and scores using those
  pre-trained parameters.

It then produces three faceted scatter plots (`CompareScores`), each with a y = x line: branch array vs branch WGS,
branch array vs `main` array, and branch array vs branch array scored with the separately pre-trained model. If the
maximum absolute difference in the pre-training comparison exceeds `max_diff_pretrain_threshold`, the workflow fails
with an error message — training the model separately should give identical scores to training it inline.

### Inputs

- **File? validationArrays**: array VCF to score with this branch. If not provided, `validationArraysMain` is used for
  both branches (i.e. only scoring changes, not imputation changes, are being tested).
- **File? validationArraysIndex**: index for `validationArrays`.
- **File validationArraysMain**: array VCF to score with the `main` branch.
- **File validationArraysIndexMain**: index for `validationArraysMain`.
- **File validationWgs**: WGS VCF for the same samples, used as the comparison truth.
- **String population_basename**: names the population-related output files.
- **File population_loadings**: PC loadings from `PerformPopulationPCA`.
- **File population_meansd**: PC means/SDs from `PerformPopulationPCA`.
- **File population_pcs**: population PCs from `PerformPopulationPCA`.
- **File pruning_sites_for_pca**: LD-pruned site list from `PerformPopulationPCA`.
- **File population_vcf**: reference population VCF (`sorted_variant_id_dataset` from `PerformPopulationPCA`).
- **Array[NamedWeightSet] named_weight_sets**: the conditions/weight sets to validate.
- **File sample_name_map**: maps sample names in the arrays to sample names in the WGS VCF, one pair per line with
  `:` as the separator.
- **String branch**: name of the branch being tested. Used for plot labels and output file names only; it does not
  affect computation.
- **Int wgs_vcf_to_plink_mem = 8**: memory (GB) for converting the WGS VCF to a plink dataset.
- **Float max_diff_pretrain_threshold = 0.0000000001**: maximum adjusted-score difference tolerated between training
  the ancestry model separately and training it as part of scoring. Exceeding it fails the workflow.

### Outputs

- **File score_comparison_branch**: `score_comparison_<branch>.png` — branch array scores vs branch WGS scores.
- **File score_comparison_main_vs_branch**: `score_comparison_main_vs_<branch>.png` — branch array scores vs `main`
  array scores.
- **File score_comparison_pretrian**: `score_comparison_pre_train_model_<branch>.png` — branch array scores vs scores
  computed with the separately pre-trained model. (Output name is misspelled in the WDL.)
- **Float pretrain_max_diff**: maximum absolute adjusted-score difference in the pre-training comparison.
- **Float branch_main_max_diff**: maximum absolute adjusted-score difference between the branch and `main`.
- **File pc_plot**: the PCA plot from the first scored condition.
- **Array[Int] n_original_weights**: per condition, the number of weights before subsetting.
- **Array[Int] n_subset_weights**: per condition, the number of weights after subsetting to the imputed array sites.
- **Array[Int] n_subset_weights_wgs**: per condition, the number of weights after further subsetting to the WGS sites.

### Docker images

`rocker/tidyverse` for the score comparison, `skwalker/imputation:with_vcftools` for `QCSites`, plus the images used
by the imported scoring and PCA tasks.

## SubsetWeightSet

Defined in [SubsetWeightSet.wdl](SubsetWeightSet.wdl). Not registered on Dockstore; it is imported by
[ValidateScoring.wdl](ValidateScoring.wdl).

### Summary

Subsets every component of a [`WeightSet`](../Structs.wdl) to a given list of site IDs, returning a new `WeightSet`.
The linear weights are filtered on their first column; the interaction weights are kept only when **both** `id_1` and
`id_2` survive; the self-exclusive sites are filtered on their `id` column, preserving the original `maxAllowed`. The
interaction weights and self-exclusive sites are only processed if they are present in the input weight set.

### Inputs

- **WeightSet weight_set**: the weight set to subset (linear weights, and optionally interaction weights and
  self-exclusive sites).
- **File sites_to_subset_to**: a file with one site ID per line; only weights whose site IDs all appear here are kept.

### Outputs

- **WeightSet subsetted_weight_set**: the subsetted weight set, with each component file prefixed `subset_`.
- **Int n_original_weights**: number of linear weights plus interaction weights before subsetting.
- **Int n_subset_weights**: number of linear weights plus interaction weights after subsetting.

**Note:** in the current `SubsetSitesBasedFile` task both `n_original.txt` and `n_subset.txt` are written from the
*subsetted* table, so `n_original_weights` and `n_subset_weights` will always be equal.

### Docker images

`rocker/tidyverse`.
