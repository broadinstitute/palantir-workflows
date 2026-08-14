# Polygenic Risk Scores (PRS)

This directory contains the WDLs, task libraries, docker definitions and python helper packages used by the DSP
Palantir team to compute, ancestry-adjust, aggregate and QC polygenic risk scores.

The core idea implemented here is:

1. **Score** a sample (or cohort) against a set of per-variant linear weights (and, optionally, pairwise interaction
   weights) using `plink2 --score`.
2. **Project** the sample onto principal components derived from a reference population (e.g. Thousand Genomes) using
   `flashPCA`.
3. **Adjust** the raw score with a model, fit on the reference population, which predicts both the mean and the
   variance of the raw score as a function of the first four PCs. The adjusted score is a z-score, and a percentile is
   derived from it with `pnorm`.
4. **Aggregate / report** the per-sample results for a lab batch into tables, plots and an HTML report.

## Contents

### Top-level workflows

| WDL | Workflow | Summary |
| --- | --- | --- |
| [ScoringPart.wdl](ScoringPart.wdl) | `ScoringImputedDataset` | The main scoring workflow: scores an imputed VCF against a weight set, projects onto reference PCs, and ancestry-adjusts the score. See [ScoringImputedDataset](#scoringimputeddataset). |
| [PerformPopulationPCA.wdl](PerformPopulationPCA.wdl) | `PerformPopulationPCA` | One-time preparation of a reference population: LD-prunes the population VCF and runs flashPCA to produce the loadings/meansd/PCs used by every downstream scoring run. See [PerformPopulationPCA](#performpopulationpca). |
| [TrainAncestryAdjustmentModel.wdl](TrainAncestryAdjustmentModel.wdl) | `TrainAncestryAdjustmentModel` | Scores the reference population with a weight set and fits the ancestry-adjustment model for that weight set. See [TrainAncestryAdjustmentModel](#trainancestryadjustmentmodel). |
| [PRSWrapper.wdl](PRSWrapper.wdl) | `PRSWrapper` | Runs `ScoringImputedDataset` for a list of conditions on one sample, applies the CKD/APOL1 adjustment where relevant, applies a reportable range and high-risk threshold, and emits a single per-sample results row. See [PRSWrapper](#prswrapper). |
| [AggregatePRSResults.wdl](AggregatePRSResults.wdl) | `AggregatePRSResults` | Aggregates `PRSWrapper` outputs for a lab batch into combined tables, a PCA plot, a score-distribution plot and an interactive HTML report. See [AggregatePRSResults](#aggregateprsresults). |
| [CKDRiskAdjustment.wdl](CKDRiskAdjustment.wdl) | `CKDRiskAdjustment` | Genotypes the APOL1 G1/G2 risk alleles and adds an APOL1 high-risk indicator to an adjusted CKD score. See [CKDRiskAdjustment](#ckdriskadjustment). |
| [PCARE.wdl](PCARE.wdl) | `PCARE` | Combined BGE (exome GVCF + imputed WGS VCF) PRS, PC projection and family history into a single linear combined risk score with low/average/high determination. See [PCARE](#pcare). |
| [PCAREAndQC.wdl](PCAREAndQC.wdl) | `PCAREAndQC` | Runs `PCARE` and then the `PRSQC` workflow from `Utilities/WDLs`. See [PCAREAndQC](#pcareandqc). |
| [ScoreBGE/ScoreBGE.wdl](ScoreBGE/ScoreBGE.wdl) | `ScoreBGE` | Scores a sample from both a WES GVCF and an imputed WGS VCF, preferring the exome genotype where it is high quality. Documented in [ScoreBGE/README.md](ScoreBGE/README.md). |
| [Validation/ValidateScoring.wdl](Validation/ValidateScoring.wdl) | `ValidateScoring` | Regression-test workflow comparing scoring on a branch against `main` and against WGS truth. Documented in [Validation/README.md](Validation/README.md). |
| [Validation/SubsetWeightSet.wdl](Validation/SubsetWeightSet.wdl) | `SubsetWeightSet` | Subsets a `WeightSet` (linear weights, interaction weights and self-exclusive sites) to a given list of site IDs. Documented in [Validation/README.md](Validation/README.md). |

### Task libraries and structs

| File | Role |
| --- | --- |
| [ScoringTasks.wdl](ScoringTasks.wdl) | Shared scoring, model-fitting and score-adjustment tasks. See [ScoringTasks](#scoringtasks). |
| [PCATasks.wdl](PCATasks.wdl) | Shared PCA / projection / plink-conversion tasks. See [PCATasks](#pcatasks). |
| [Structs.wdl](Structs.wdl) | Shared WDL structs. See [Structs](#structs). |

### How the workflows chain together

```
                        PerformPopulationPCA          (once per reference population + array platform)
                                 |
        population_loadings, population_meansd, population_pcs, pruning_sites_for_pca,
        sorted_variant_id_dataset (population_vcf)
                                 |
                                 v
                 TrainAncestryAdjustmentModel          (once per weight set; also callable inline)
                                 |
                    fitted_model_params + sites_used_in_scoring
                          (= AncestryAdjustmentModelParams)
                                 |
                                 v
     imputed VCF ----> ScoringImputedDataset (ScoringPart.wdl) ----> adjusted_array_scores, pc_projection
                                 ^
                                 |
                            PRSWrapper           (scatter over conditions, one sample)
                                 |
                    results.csv, pcs, missing_sites_shifts
                                 |
                                 v
                         AggregatePRSResults      (one lab batch)
                                 |
              batch tables, PCA plot, score distribution, HTML report
```

`ScoringImputedDataset` can obtain the ancestry-adjustment model in two mutually exclusive ways: either you pass a
pre-fitted `fitted_model_params_and_sites`, or you pass a `population_vcf` and the workflow calls
`TrainAncestryAdjustmentModel` itself. Production use (`PRSWrapper`) passes pre-fitted parameters.

`PCARE` / `PCAREAndQC` are a separate, simpler track: they do not use the ancestry-adjustment model at all. They
combine a BGE PRS with PCs and family history into a single linear risk score.

## Structs

[Structs.wdl](Structs.wdl) defines the structs shared across these workflows. `Structs.wdl` is imported (often with a
`#!UnusedImport` womtool-lint suppression) by every workflow that takes one of these as an input.

- **ReferencePanelContig**: `File vcf`, `File vcf_index`, `File bcf`, `File bcf_index`, `File m3vcf`, `String contig`.
  A per-contig imputation reference panel bundle. Not used by any workflow in this directory; retained for
  imputation workflows that import these structs.
- **AncestryAdjustmentModelParams**: `File fitted_model_params` (the 10 fitted mean/variance coefficients written by
  `ScoringTasks.TrainAncestryModel`) and `File sites_used_in_scoring` (the sites the model was trained over). Passing
  both together guarantees that a score being adjusted was computed over the same site set the model was fit on.
- **SelfExclusiveSites**: `File sites` (must have columns `id`, `chrom`, `pos`, and is read with an `allele` column by
  `AddInteractionTermsToScore`) and `Int maxAllowed`. Used to suppress interaction terms in samples that carry more
  than `maxAllowed` of the listed effect alleles.
- **WeightSet**: `File linear_weights` (standard PRS weights file), optional `File interaction_weights` (must have
  columns `id_1`, `id_2`, `chrom_1`, `chrom_2`, `pos_1`, `pos_2`, `allele_1`, `allele_2`, `weight`; order not
  important) and optional `SelfExclusiveSites interaction_self_exclusive_sites`.
- **NamedWeightSet**: `String condition_name` + `WeightSet weight_set`. The condition name is used to name outputs and
  to label columns in the aggregated results.
- **PRSWrapperConditionResource**: `Boolean score_condition`, `Float percentile_threshold`, `NamedWeightSet
  named_weight_set`, `AncestryAdjustmentModelParams ancestry_model_params_and_sites`. One element per condition in the
  `PRSWrapper` input array.

## ScoringImputedDataset

Defined in [ScoringPart.wdl](ScoringPart.wdl); registered on Dockstore as **PRScoringWorkflow**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/ScoringPart.html) · [open locally](../docs/viz/PRS/ScoringPart.html)

### Summary

The main user-facing scoring workflow. Given an imputed VCF and a `NamedWeightSet`, it:

- determines the chromosome encoding to hand to plink from the weights file (`DetermineChromosomeEncoding`);
- scores the VCF, either with `plink2 --score` (`ScoringTasks.ScoreVcf`) or, if `use_bge_scoring` is set, with the
  [ScoreBGE](ScoreBGE/README.md) workflow against a WES GVCF plus the imputed VCF;
- optionally adds pairwise interaction terms if the weight set has `interaction_weights`;
- if `adjustScores` (default), converts the VCF to a plink dataset restricted to the PCA pruning sites, projects it
  onto the population PCs with flashPCA, and adjusts the raw score with the fitted model;
- compares the sites actually scored against the sites used when the model was trained. If sites are missing, it
  computes the maximum possible upward and downward shift in the raw score (2x the weight for each missing linear
  site) and reports the resulting range of adjusted scores/percentiles.

Note that variant IDs in the weights file, the imputed VCF and the population files must match exactly. IDs are set by
plink from `--set-all-var-ids`, using `chr:pos:ref:alt` when `use_ref_alt_for_ids` is true and `chr:pos:allele1:allele2`
(sorted alleles) otherwise.

The workflow uses `ErrorWithMessage` tasks to fail fast on invalid input combinations: `adjustScores` requires all of
`population_loadings`, `population_meansd`, `population_pcs` and `pruning_sites_for_pca`, and requires exactly one of
`fitted_model_params_and_sites` or `population_vcf`; `use_bge_scoring` requires all of `bge_wes_gvcf`,
`bge_wes_gvcf_index` and `ref_dict`.

### Inputs

- **NamedWeightSet named_weight_set**: the condition name and weight set (linear weights, and optionally interaction
  weights and self-exclusive sites) to score with.
- **File imputed_array_vcf**: the imputed VCF to score, and optionally to project for PCA. Variant IDs must match
  those in the weights file.
- **File? imputed_array_vcf_index**: index for `imputed_array_vcf`. Required when `use_bge_scoring` is true.
- **Boolean use_bge_scoring = false**: score with the `ScoreBGE` workflow (WES GVCF + imputed WGS VCF) instead of
  `plink2 --score` on the VCF alone.
- **File? bge_wes_gvcf**: WES GVCF for BGE scoring. Required if `use_bge_scoring`.
- **File? bge_wes_gvcf_index**: index for `bge_wes_gvcf`. Required if `use_bge_scoring`.
- **File? ref_dict**: reference sequence dictionary, used by `ScoreBGE` to order weights by contig. Required if
  `use_bge_scoring`.
- **Int scoring_mem = 16**: base memory (GB) for the `ScoreVcf` task; plink is given 75% of it.
- **Int vcf_to_plink_mem = 8**: memory (GB) for `ExtractIDsPlink` and `ArrayVcfToPlinkDataset`.
- **String? population_basename**: used to name the population scoring outputs when the ancestry model is trained
  inline. Required if `population_vcf` is provided.
- **String basename**: names the array scoring outputs and the array projection files.
- **File? population_loadings**: PC loadings from `PerformPopulationPCA`. Required if `adjustScores`.
- **File? population_meansd**: PC means/SDs from `PerformPopulationPCA`. Required if `adjustScores`.
- **File? population_pcs**: population PCs from `PerformPopulationPCA`. Required if `adjustScores`.
- **File? pruning_sites_for_pca**: LD-pruned site list from `PerformPopulationPCA`. Required if `adjustScores`.
- **File? population_vcf**: reference population VCF (the `sorted_variant_id_dataset` output of
  `PerformPopulationPCA`). If given, the ancestry-adjustment model is trained inline. Mutually exclusive with
  `fitted_model_params_and_sites`.
- **AncestryAdjustmentModelParams? fitted_model_params_and_sites**: a previously fitted model and the sites it was
  trained over. Mutually exclusive with `population_vcf`.
- **String? columns_for_scoring**: passed through to `plink2 --score` as extra arguments. plink expects the first three
  columns of the weights file to be variant ID, effect allele and effect weight; if that is not the case, give the
  column numbers in that order, e.g. `"11 12 13"`.
- **Boolean redoPCA = false**: when `population_vcf` is provided, re-run `PerformPCA` on the population restricted to
  the sites present in the target VCF instead of using the supplied loadings/meansd/PCs.
- **Boolean adjustScores = true**: perform PC projection and ancestry adjustment. If false, only raw scores are
  produced.
- **Boolean use_ref_alt_for_ids = false**: build variant IDs as `chr:pos:ref:alt` rather than `chr:pos:allele1:allele2`
  with sorted alleles.

### Outputs

- **File? pc_projection**: the target sample(s) projected onto the population PCs (`projections.txt` from flashPCA).
  Present when `adjustScores`.
- **File raw_scores**: the unadjusted score table; the interaction-augmented table if the weight set has interaction
  weights, otherwise the plink `.sscore` (or `ScoreBGE` `.score`) table.
- **File? pc_plot**: PNG scatter of PC1 vs PC2 with the reference population and the target sample(s) overlaid.
- **File? adjusted_population_scores**: adjusted scores for the reference population; only produced when the model was
  trained inline from `population_vcf`.
- **File? adjusted_array_scores**: the adjusted scores table, with `adjusted_score` (z-score) and `percentile` columns.
- **Boolean? fit_converged**: whether the `optim` BFGS fit of the ancestry model converged; only produced when the
  model was trained inline.
- **Int? n_missing_sites_from_training**: number of sites used in model training that were not scored in this run.
- **File? missing_sites_shifted_scores**: per-sample table with `condition`, `n_missing_sites`, `adjusted_score`,
  `percentile` and the potential high/low adjusted score and percentile implied by the missing sites.

## PerformPopulationPCA

Defined in [PerformPopulationPCA.wdl](PerformPopulationPCA.wdl); registered on Dockstore as **PerformPopulationPCA**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/PerformPopulationPCA.html) · [open locally](../docs/viz/PRS/PerformPopulationPCA.html)

### Summary

Run this once whenever you adopt a new reference population dataset (e.g. Thousand Genomes) or a new array platform.
It normalises the population VCF's variant IDs, restricts it to the sites present in the imputed arrays, LD-prunes,
and runs flashPCA. The outputs are the inputs that `ScoringImputedDataset` needs for PC projection and adjustment.

Steps: split multiallelics and set IDs to `CHROM:POS:REF:FIRST_ALT` (`SeparateMultiallelics`, bcftools); sort the
alleles within each ID so they are order-independent (`SortVariantIds`); do the same for each imputed array VCF and
select only `TYPED || TYPED_ONLY` sites from it (`UpdateVariantIds`, `SelectTypedSites`, GATK `SelectVariants`);
intersect the population VCF down to those intervals (`SubsetToArrayVCF`); LD-prune with plink2
(`--geno 0.05 --hwe 1e-10 --maf 0.01 --indep-pairwise 1000 50 0.2`, excluding chromosome X, intersected with the
original-array snplists) (`LDPruning`); run `PCATasks.PerformPCA`; and finally run flashPCA's `--check` on the result
(`CheckPCA`) so you can inspect the mean squared error in the task log (flashPCA documentation suggests looking for
`<1e-8`).

`LDPruneToSites` is also defined in this file but is not called by the workflow; use it if you already have a list of
LD-pruned sites to prune to.

### Inputs

- **File population_vcf**: the reference population VCF (e.g. Thousand Genomes).
- **File population_vcf_index**: index for `population_vcf`.
- **String basename**: prefix for the output file names.
- **Array[File] imputed_array_vcfs**: imputed array VCFs. These are limited to TYPED and TYPED_ONLY sites before LD
  pruning, and also constrain the population to sites usable for scoring correction.
- **Array[File] original_array_vcfs**: the pre-imputation array VCFs; sites are selected from these with plink2
  (`--geno 0.001 --snps-only`) and intersected into the pruning step.
- **Array[File]? subset_to_sites**: optional additional site lists to intersect with during LD pruning.
- **String? chromosome_encoding**: plink2 `--output-chr` encoding (see the
  [plink2 docs](https://www.cog-genomics.org/plink/2.0/data#irreg_output)).

### Outputs

- **File population_loadings**: flashPCA PC loadings (`<basename>.pc.loadings`).
- **File population_meansd**: flashPCA per-variant means and SDs (`<basename>.pc.meansd`).
- **File population_pcs**: the population PCs (`<basename>.pc`).
- **File pruning_sites_for_pca**: the LD-pruned site list (`<basename>.prune.in`).
- **File sorted_variant_id_dataset**: the population VCF with sorted, normalised variant IDs. **This** is what you
  should pass as `population_vcf` to `ScoringImputedDataset` / `TrainAncestryAdjustmentModel`, so that IDs match.
- **File sorted_variant_id_dataset_index**: index for `sorted_variant_id_dataset`.

### Docker images

`biocontainers/bcftools:v1.9-1-deb_cv1`, `us.gcr.io/broad-gatk/gatk:4.1.9.0`, `skwalker/plink2:first`,
`skwalker/imputation:with_vcftools`, `skwalker/flashpca:v1`, and (via `PCATasks.PerformPCA`)
`us.gcr.io/broad-dsde-methods/flashpca_docker`.

## TrainAncestryAdjustmentModel

Defined in [TrainAncestryAdjustmentModel.wdl](TrainAncestryAdjustmentModel.wdl); registered on Dockstore as
**TrainAncestryAdjustmentModel**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/TrainAncestryAdjustmentModel.html) · [open locally](../docs/viz/PRS/TrainAncestryAdjustmentModel.html)

### Summary

Scores the reference population with a weight set and fits the ancestry-adjustment model for it. The fit is done in
`ScoringTasks.TrainAncestryModel`: a Gaussian GLM of `SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4` and a Gamma(log-link) GLM of
the squared residuals give starting values, and then the full 10-parameter mean-and-variance likelihood is maximised
with `optim(method = "BFGS")`. Adjusted scores are `(SCORE1_SUM - mu(PCs)) / sqrt(sigma2(PCs))` and percentiles are
`pnorm` of that.

Called on its own to produce reusable `AncestryAdjustmentModelParams`, or called inline by `ScoringImputedDataset`
when a `population_vcf` is provided.

### Inputs

- **NamedWeightSet named_weight_set**: condition name and weight set to fit the model for.
- **File population_pcs**: population PCs, from `PerformPopulationPCA`.
- **File population_vcf**: the reference population VCF, from `PerformPopulationPCA`
  (`sorted_variant_id_dataset`). Variant IDs must match those in the weights file.
- **String population_basename**: prefix for the population scoring and fitted-model output files.
- **File? sites**: restrict scoring to this set of site IDs (typically the sites present in the target VCF).

### Outputs

- **File fitted_params**: TSV of the 10 fitted model parameters (`<condition>_<population_basename>_fitted_model_params.tsv`).
- **File sites_used_in_scoring**: the sites the model was trained over (linear scoring sites, unioned with interaction
  sites if interaction weights were supplied).
- **File adjusted_population_scores**: the reference population's adjusted scores and percentiles.
- **Boolean fit_converged**: whether the BFGS optimisation converged.

## PRSWrapper

Defined in [PRSWrapper.wdl](PRSWrapper.wdl); registered on Dockstore as **PRSWrapper**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/PRSWrapper.html) · [open locally](../docs/viz/PRS/PRSWrapper.html)

### Summary

The production per-sample entry point. It scatters over a list of `PRSWrapperConditionResource`s (one per condition),
calling `ScoringImputedDataset` for each condition whose `score_condition` is true, and produces a single-row CSV of
results for the sample. For the condition literally named `ckd`, [CKDRiskAdjustment](#ckdriskadjustment) is applied to
the adjusted score before reporting.

For each scored condition it emits `<condition>_raw`, `<condition>_adjusted`, `<condition>_percentile`,
`<condition>_risk` and `<condition>_reason_not_resulted` columns. `risk` is `HIGH` when the percentile exceeds the
condition's `percentile_threshold`, otherwise `NOT_HIGH`. If `|adjusted_score| > z_score_reportable_range`, all four
value columns become `NOT_RESULTED` and `reason_not_resulted` records whether the z-score was above or below the
range. Conditions with `score_condition = false` get a row of `NA`s from `CreateUnscoredResult`.

`SelectValuesOfInterest` and `CheckZScoreAgainstReportableRange` both require the score table to contain exactly one
row, i.e. this workflow is for a single sample at a time. `SelectValuesOfInterest` also fails if the `IID` in the
score table does not match the provided `sample_id`.

### Inputs

- **Array[PRSWrapperConditionResource] condition_resources**: one entry per condition, giving whether to score it, the
  high-risk percentile threshold, the named weight set and the pre-fitted ancestry model parameters.
- **File ckd_risk_alleles**: APOL1 G1/G2 risk allele site list, passed to `CKDRiskAdjustment`.
- **Float z_score_reportable_range**: absolute adjusted-score (z-score) limit beyond which results are reported as
  `NOT_RESULTED`.
- **File vcf**: the imputed VCF for the sample.
- **Boolean use_bge_scoring = false**: use BGE (WES GVCF + imputed VCF) scoring in `ScoringImputedDataset`.
- **File? bge_wes_gvcf**: WES GVCF, required when `use_bge_scoring`.
- **File? bge_wes_gvcf_index**: index for `bge_wes_gvcf`.
- **File? ref_dict**: reference sequence dictionary, required when `use_bge_scoring`.
- **String sample_id**: sample ID; must match the `IID` in the score tables.
- **String lab_batch_id**: lab batch identifier, added as a column and passed through as an output.
- **Boolean is_control_sample_in**: whether this sample is the batch's control sample; added as a column and passed
  through as an output.
- **Boolean redoPCA = false**: forwarded to `ScoringImputedDataset`.
- **File population_loadings**: PC loadings from `PerformPopulationPCA`.
- **File population_meansd**: PC means/SDs from `PerformPopulationPCA`.
- **File population_pcs**: population PCs from `PerformPopulationPCA`.
- **File pruning_sites_for_pca**: LD-pruned site list from `PerformPopulationPCA`.
- **Boolean use_ref_alt_for_ids = false**: forwarded to `ScoringImputedDataset` and `CKDRiskAdjustment`.
- **Int? vcf_to_plink_mem**: forwarded to `ScoringImputedDataset` (`vcf_to_plink_mem`) and `CKDRiskAdjustment`
  (`mem_plink`).

### Outputs

- **File results**: single-row CSV with `sample_id`, `lab_batch`, `is_control_sample`, and the five columns per
  condition described above.
- **File pcs**: the PC projection for this sample (taken from the first scored condition).
- **String lab_batch**: passthrough of `lab_batch_id`.
- **Boolean is_control_sample**: passthrough of `is_control_sample_in`.
- **File missing_sites_shifts**: `<lab_batch>.missing_sites_shifts.tsv`, the concatenation of the per-condition
  missing-site shift tables.

### Docker images

`rocker/tidyverse` (pinned by digest) for all wrapper tasks; the scoring tasks use their own images (see
[ScoringTasks](#scoringtasks)).

## AggregatePRSResults

Defined in [AggregatePRSResults.wdl](AggregatePRSResults.wdl); registered on Dockstore as **AggregatePRSResults**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/AggregatePRSResults.html) · [open locally](../docs/viz/PRS/AggregatePRSResults.html)

### Summary

Aggregates the per-sample `PRSWrapper` outputs for one lab batch. It joins the results with the PC projections,
validates that all inputs come from a single lab batch matching `lab_batch` and that exactly one control sample is
present, writes long/wide/summary tables, plots the batch score distribution against a standard normal, plots the
batch PCs over the reference population PCs, and renders an interactive HTML report with `rmarkdown`.

The HTML report contains: a control-sample table comparing the observed control against `expected_control_results`
(deltas coloured red when missing or larger than `control_sample_diff_threshold`); a batch summary; a table of the
condition groups actually scored per sample checked against `allowed_condition_groups`; a table of samples high risk
for multiple conditions with a significance computed from the per-condition thresholds under an
all-conditions-uncorrelated assumption; an interactive plotly score-distribution plot; an interactive plotly PCA plot;
a DataTables view of individual sample results; and a missing-sites table.

Outputs are prefixed `<lab_batch>` when `group_n` is 1 and `<lab_batch>_group_<group_n>` otherwise.

### Inputs

- **Array[File] results**: the `results` CSVs from `PRSWrapper`, one per sample.
- **Array[File] target_pc_projections**: the `pcs` files from `PRSWrapper`, one per sample.
- **Array[File] missing_sites_shifts**: the `missing_sites_shifts` files from `PRSWrapper`, one per sample.
- **File high_risk_thresholds**: TSV with `condition` and `threshold` columns; used to compute multi-high-risk
  significance in the report.
- **File population_pc_projections**: reference population PCs, for the PCA plots.
- **String population_name = "Reference Population"**: label for the reference population in the PCA plots.
- **File expected_control_results**: CSV of the expected control sample values, compared against the observed control
  in the report.
- **File allowed_condition_groups**: TSV with `group` and `condition` columns defining which sets of conditions are
  allowed to be scored together; samples whose scored set is not in this list are flagged red in the report.
- **String lab_batch**: expected lab batch ID; the workflow fails if the inputs disagree.
- **Int group_n**: group number within the lab batch; values > 1 add a `_group_<n>` suffix to output names and to the
  report title.
- **Float control_sample_diff_threshold**: absolute difference above which a control-sample delta is coloured red.

### Outputs

- **File batch_all_results**: all per-sample results joined with PCs.
- **File batch_control_results**: the control sample's row.
- **File batch_summarised_results**: per-condition means, sample counts and high / not-high / not-resulted counts.
- **File batch_missing_sites_shifts**: concatenated missing-site shift tables.
- **File score_distribution**: PNG density plot of adjusted z-scores per condition against a standard normal.
- **File pc_plot**: PNG PC1 vs PC2 plot of the batch over the reference population.
- **File report**: the rendered HTML report.
- **File batch_pcs**: the combined batch PC table.

### Docker images

`rocker/tidyverse` (pinned by digest) for aggregation and plotting;
`us.gcr.io/broad-dsde-methods/tidyverse_kableextra_docker` for the report (it adds `kableExtra`, `plotly` and `DT`).

## CKDRiskAdjustment

Defined in [CKDRiskAdjustment.wdl](CKDRiskAdjustment.wdl); registered on Dockstore as **CKDRiskAdjustment**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/CKDRiskAdjustment.html) · [open locally](../docs/viz/PRS/CKDRiskAdjustment.html)

### Summary

Applies an APOL1 adjustment on top of an already-adjusted chronic kidney disease score. It genotypes the G1/G2 risk
alleles out of the VCF with `plink2 --export A --export-allele`, calls a sample APOL1 high risk when the sum of the
larger of the two G1 allele counts and the G2 allele count is at least 2, and then adds 1 to the adjusted score for
high-risk samples, recomputing the percentile with `pnorm`.

### Inputs

- **File adjustedScores**: the adjusted scores table from `ScoringImputedDataset`.
- **File vcf**: the imputed VCF for the sample.
- **File risk_alleles**: site list for the APOL1 G1 (two sites) and G2 risk alleles, in a plink `--export-allele`
  compatible format. Also used to determine the chromosome encoding.
- **Boolean use_ref_alt_for_ids = false**: build variant IDs as `chr:pos:ref:alt` rather than sorted alleles.
- **Int mem_plink = 8**: memory (GB) for the plink task; plink is given 75% of it.

### Outputs

- **File adjusted_scores_with_apol1**: the adjusted scores table with an `apol1_high_risk` column and the
  APOL1-adjusted `adjusted_score` and `percentile`.

### Docker images

`us.gcr.io/broad-dsde-methods/plink2_docker` and `rocker/tidyverse:4.1.0`.

## PCARE

Defined in [PCARE.wdl](PCARE.wdl); registered on Dockstore as **PCARE**. Covered by a WATT test (see [Testing](#testing)).

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/PCARE.html) · [open locally](../docs/viz/PRS/PCARE.html)

### Summary

Computes a combined risk value from a BGE PRS, two principal components and a family history indicator. It runs
[ScoreBGE](ScoreBGE/README.md) on the exome GVCF plus the imputed WGS VCF, converts the imputed WGS VCF to a plink
dataset restricted to `pc_sites` and projects it onto the supplied PC loadings (with flashPCA `--div none`), then
computes

```
combined_risk_score = prs_beta*SCORE1_SUM + fam_hist_beta*fam_hist + pc1_beta*PC1 + pc2_beta*PC2
```

and classifies it as `low`, `average` or `high` using the two thresholds. All of the betas and thresholds are inputs;
the defaults for the thresholds (19.69 / 20.38) are the only values baked into the WDL.

### Inputs

- **File imputed_wgs_vcf**: imputed WGS VCF; used both for scoring and for PC projection.
- **File imputed_wgs_vcf_index**: index for `imputed_wgs_vcf`.
- **File exome_gvcf**: WES GVCF for `ScoreBGE`.
- **File exome_gvcf_index**: index for `exome_gvcf`.
- **File prs_weights**: weights file in the `ScoreBGE` format (`contig`, `position`, `ref`, `alt`, `effect_allele`,
  `weight`).
- **Array[String]? sample_names**: restrict scoring to these samples; if absent, all samples are scored.
- **File fam_history**: TSV with a `sample_id` column and a `fam_hist` column.
- **String basename**: prefix for output file names.
- **File pc_loadings**: PC loadings for projection.
- **File pc_meansd**: PC means/SDs for projection.
- **File pc_sites**: sites to restrict the plink dataset to before projection.
- **Float prs_beta**: coefficient on the PRS score in the combined risk score.
- **Float fam_hist_beta**: coefficient on the family history indicator.
- **Float pc1_beta**: coefficient on PC1.
- **Float pc2_beta**: coefficient on PC2.
- **Float risk_determination_threshold_low_average = 19.69**: below this the risk determination is `low`.
- **Float risk_determination_threshold_average_high = 20.38**: above this the risk determination is `high`.
- **Boolean use_ref_alt_for_ids = true**: build variant IDs as `chr:pos:ref:alt`.
- **String chromosome_encoding = "chrMT"**: plink2 `--output-chr` encoding.
- **Int mem_gb_array_vcf_to_plink = 16**: memory (GB) for the VCF-to-plink conversion.
- **File ref_dict**: reference sequence dictionary, used by `ScoreBGE`.

### Outputs

- **File full_risk**: `<basename>_full_risk.tsv`, with columns `sample_id`, `prs_score`, `pc1`, `pc2`,
  `family_history`, `combined_risk_score` and `risk_determination`.

### Docker images

`us.gcr.io/broad-dsde-methods/plink2_docker` and `us.gcr.io/broad-dsde-methods/flashpca_docker` (via `PCATasks`),
`us.gcr.io/broad-dsde-methods/python-data-slim:1.0` for the risk computation, and the ScoreBGE image.

## PCAREAndQC

Defined in [PCAREAndQC.wdl](PCAREAndQC.wdl); registered on Dockstore as **PCAREAndQC**.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/PRS/PCAREAndQC.html) · [open locally](../docs/viz/PRS/PCAREAndQC.html)

### Summary

Runs [PCARE](#pcare) and then the `PRSQC` workflow from
[../Utilities/WDLs/PRSQC.wdl](../Utilities/WDLs/PRSQC.wdl), which checks the score against an acceptable range and
checks whether the sample's (PC1, PC2) falls inside an alphashape bounding polygon fit to a training population.

### Inputs

All of the [PCARE](#pcare) inputs, plus:

- **File acceptable_range**: score range definition used by `PRSQC.CheckScores`.
- **File alphashape**: bounding polygon over (PC1, PC2) from a training set, generated with the
  [alphashape](https://github.com/bellockk/alphashape) library. See the comments in
  [../Utilities/WDLs/PRSQC.wdl](../Utilities/WDLs/PRSQC.wdl) for details.
- **Float distance_threshold**: how far outside the alphashape a point may fall and still pass.

### Outputs

- **File full_risk**: passthrough of the `PCARE` output.
- **Boolean qc_passed**: true when both the score range check and the PCA shape check pass.
- **Boolean pcs_within_shape**: whether the sample's PCs fall within the alphashape (plus `distance_threshold`).
- **File pca_qc_plot**: plot of the sample's PCs against the alphashape.

## ScoringTasks

[ScoringTasks.wdl](ScoringTasks.wdl) is a task library, not a workflow. It holds every task involved in computing,
shifting, fitting and adjusting scores, and is imported by `ScoringPart.wdl`, `TrainAncestryAdjustmentModel.wdl`,
`CKDRiskAdjustment.wdl` and `Validation/ValidateScoring.wdl`.

- **ScoreVcf**: runs `plink2 --score` over a VCF with a weights file, using `DS` dosages by default
  (`use_dosage_annotation`), with `no-mean-imputation` and `ignore-dup-ids`; emits the `.sscore` table, the plink log
  and the list of sites scored. Docker: `us.gcr.io/broad-dsde-methods/plink2_docker`.
- **AddInteractionTermsToScore**: python/cyvcf2 task that counts effect alleles at each interaction site, adds the
  product of the paired allele counts times the interaction weight to `SCORE1_SUM`, and honours
  `SelfExclusiveSites` (samples carrying more than `maxAllowed` of the listed alleles get no interaction terms).
  Docker: `us.gcr.io/broad-dsde-methods/imputation_interaction_python`.
- **CheckWeightsCoverSitesUsedInTraining**: fails the workflow if any site used to train the model is absent from the
  weight set. Docker: `python:3.9.10`.
- **CompareScoredSitesToSitesUsedInTraining**: computes the sites used in training but not scored, and the maximum
  possible upward and downward shift to the raw score they could account for. Docker: `python:3.9.10`.
- **CombineScoringSites**: unions the linear and interaction scoring site lists. Docker: `ubuntu:20.04`.
- **AddShiftToRawScores**: adds a fixed shift to `SCORE1_SUM`. Docker: `rocker/tidyverse:4.1.0`.
- **CombineMissingSitesAdjustedScores**: joins the adjusted, shifted-up and shifted-down score tables into the
  `missing_sites_shifted_scores` report. Docker: `rocker/tidyverse:4.1.0`.
- **TrainAncestryModel**: fits the mean-and-variance ancestry-adjustment model on the reference population (see
  [TrainAncestryAdjustmentModel](#trainancestryadjustmentmodel)). Docker: `rocker/tidyverse` (pinned by digest).
- **AdjustScores**: applies a fitted model to a score table and a PC table, producing `adjusted_score` and
  `percentile`. Docker: `rocker/tidyverse` (pinned by digest).
- **MakePCAPlot**: PC1 vs PC2 ggplot of the population and target PCs. Docker: `rocker/tidyverse` (pinned by digest).
- **ExtractIDsPlink**: writes the variant ID list of a VCF using plink's ID conventions, excluding duplicates.
  Docker: `us.gcr.io/broad-dsde-methods/plink2_docker`.
- **DetermineChromosomeEncoding**: inspects the weights file to decide which plink `--output-chr` encoding
  (`MT` / `M` / `chrM` / `chrMT`) to use. Docker: `python:3.9.10`.

## PCATasks

[PCATasks.wdl](PCATasks.wdl) is a task library holding the PCA and plink-conversion tasks shared by
`PerformPopulationPCA.wdl`, `ScoringPart.wdl`, `PCARE.wdl` and `Validation/ValidateScoring.wdl`.

- **PerformPCA**: runs `flashpca -d 20` on a plink bed/bim/fam dataset, emitting PCs, variance explained, loadings,
  means/SDs, eigenvectors and eigenvalues. Docker: `us.gcr.io/broad-dsde-methods/flashpca_docker`.
- **ProjectArray**: projects a plink dataset onto existing PC loadings/means with `flashpca --project`. It first
  checks that the `.bim` IDs, the loadings IDs and the meansd IDs are identical and fails with an explanatory message
  otherwise. Accepts an optional `divisor` passed to flashPCA's `--div`. Docker:
  `us.gcr.io/broad-dsde-methods/flashpca_docker`.
- **ArrayVcfToPlinkDataset**: converts a VCF to plink bed/bim/fam, intersecting with the pruning sites (and optionally
  a second site list), setting variant IDs and de-duplicating. Docker: `us.gcr.io/broad-dsde-methods/plink2_docker`.

## Docker images

Most tasks use either public images (`rocker/tidyverse`, `python:3.9.10`, `ubuntu:20.04`, `biocontainers/bcftools`,
`us.gcr.io/broad-gatk/gatk`) or images built from the Dockerfiles in this directory and pushed to
`us.gcr.io/broad-dsde-methods`. WDL tasks generally pin these by digest.

| Directory | Image | Contents |
| --- | --- | --- |
| [flashpca_docker](flashpca_docker/Dockerfile) | `us.gcr.io/broad-dsde-methods/flashpca_docker` | Ubuntu bionic with R, Eigen, Boost and Spectra 0.8.1, building [github.com/kachulis/flashpca](https://github.com/kachulis/flashpca) from the `ck_project_single_sample` branch. Used by `PCATasks.PerformPCA` and `PCATasks.ProjectArray`. |
| [plink2_docker](plink2_docker/Dockerfile) | `us.gcr.io/broad-dsde-methods/plink2_docker` | plink2 (`PLINK_VERSION` build arg, default `v2.00a2.3`) built from source on Ubuntu (`UBUNTU_VERSION` build arg, default 24.04), with the binary at `/plink2`. Used by every plink task. |
| [tidyverse_kableextra_docker](tidyverse_kableextra_docker/Dockerfile) | `us.gcr.io/broad-dsde-methods/tidyverse_kableextra_docker` | `rocker/tidyverse` plus `kableExtra`, `plotly` and `DT`. Used by `AggregatePRSResults.BuildHTMLReport`. Has its own wrapper script, [build_push_r_kableextra_docker.sh](tidyverse_kableextra_docker/build_push_r_kableextra_docker.sh), which records the image version and calls `build_push_docker.sh`. |
| [imputation_interaction_python](imputation_interaction_python/Dockerfile) | `us.gcr.io/broad-dsde-methods/imputation_interaction_python` | Ubuntu with python3, tabix, `cyvcf2` and `pandas`. Used by `ScoringTasks.AddInteractionTermsToScore`. |
| [minimac3_docker](minimac3_docker/Dockerfile) | `us.gcr.io/broad-dsde-methods/minimac3_docker` | Ubuntu with the Minimac3 executable downloaded from the University of Michigan share. Not referenced by any WDL in this directory. **Note:** this directory is present in the working tree but is not tracked in git. |
| [ScoreBGE](ScoreBGE/Dockerfile) | `us.gcr.io/broad-dsde-methods/palantir-workflows-score-bge` | `python-data-slim-pysam` with `ScoreBGE.py` added at `/ScoreBGE.py`. See [ScoreBGE/README.md](ScoreBGE/README.md). |

### build_push_docker.sh

[build_push_docker.sh](build_push_docker.sh) builds and pushes an image from any of the Dockerfile directories in the
repo. It can be called from anywhere in the repo: it finds the named directory with `find ..`, tags the image as
`us.gcr.io/broad-dsde-methods/<directory-name>:<tag>`, prompts for confirmation, then builds and (unless
`--no-push`) pushes.

```
./build_push_docker.sh --directory <dockerfile-dir> --image-version-tag <tag> [options]
```

Options: `--directory/-d` (required), `--image-version-tag/-i` (required), `--ubuntu-version/-u` (passed as the
`UBUNTU_VERSION` build arg, default `20.04`), `--no-push/-p`, `--dry-run/-r`, `--no-cache/-c`, `--help/-h`.

Note that the script's default Ubuntu version (20.04) differs from the defaults declared inside
`plink2_docker/Dockerfile` and `imputation_interaction_python/Dockerfile` (24.04); the script always passes the build
arg, so the script's value wins unless you override it with `-u`.

## Python packages

Three python packages live alongside the WDLs. They are Terra/FireCloud helpers rather than parts of any WDL, and are
installed with `pip install` from their own directories (each has a `pyproject.toml`).

- **[CreateAggregationSets](CreateAggregationSets/CreateAggregationSets.py)**: builds Terra sample sets for
  `AggregatePRSResults` from the samples in a workspace, grouping by `lab_batch`, tracking which sets have already
  been delivered and creating new `_group_N` sets for reworked or late-arriving samples. Run as
  `CreateAggregationSets --workspace_namespace <ns> --workspace_name <name>`. Dependencies: `firecloud`, `pytz`.
  Tests in [CreateAggregationSets/tests](CreateAggregationSets/tests).
- **[ManualQCPRS](ManualQCPRS/ManualQCPRS.py)**: an `ipywidgets` GUI, intended to be run in a Terra notebook, for
  manually reviewing PRS batch results, editing/invalidating results and confirming delivery. Dependencies:
  `firecloud`, `pytz`, `ipywidgets`, `pandas`, `fsspec`, `gcsfs`. Tests in [ManualQCPRS/tests](ManualQCPRS/tests).
- **[Liftover/LiftoverSites](Liftover/LiftoverSites/LiftoverSites/LiftoverSites.py)**: lifts a weights or sites file
  over to another reference using Hail plus GATK/Picard `LiftoverVcf`, subsets the result to a reference panel sites
  VCF, and writes it out in either the eMERGE weights format (`CHR:BP:REF:ALT`, `effect_allele`, `weight`) or the
  `ScoreBGE` format. Dependencies: `hail`, `pandas`.

None of these three has its own README; the descriptions above are derived from the source.

## Dockstore registration

The following workflows in this directory are registered in [../.dockstore.yml](../.dockstore.yml):

| Dockstore name | Descriptor |
| --- | --- |
| PRScoringWorkflow | `/PRS/ScoringPart.wdl` |
| PerformPopulationPCA | `/PRS/PerformPopulationPCA.wdl` |
| TrainAncestryAdjustmentModel | `/PRS/TrainAncestryAdjustmentModel.wdl` |
| PRSWrapper | `/PRS/PRSWrapper.wdl` |
| AggregatePRSResults | `/PRS/AggregatePRSResults.wdl` |
| CKDRiskAdjustment | `/PRS/CKDRiskAdjustment.wdl` |
| PCARE | `/PRS/PCARE.wdl` |
| PCAREAndQC | `/PRS/PCAREAndQC.wdl` |
| ScoreBGE | `/PRS/ScoreBGE/ScoreBGE.wdl` |
| ValidateScoring | `/PRS/Validation/ValidateScoring.wdl` |

`SubsetWeightSet.wdl` is not registered on its own; it is imported by `ValidateScoring.wdl`.

## Testing

Two workflows in this directory have automated tests configured in [../test/watt_config.yml](../test/watt_config.yml):

- `AggregatePRSResults`, with inputs `/test/AggregatePRSResults/test_inputs.json` and expected outputs
  `/test/AggregatePRSResults/test_outputs.json`.
- `PCARE`, with inputs `/test/PCARE/test_inputs.json` and expected outputs `/test/PCARE/test_outputs.json`.

In addition, [Validation/ValidateScoring.wdl](Validation/ValidateScoring.wdl) is a manually-run regression workflow
that compares scoring on a development branch to scoring on `main` and to WGS-derived scores; see
[Validation/README.md](Validation/README.md). The python packages have their own unit tests under their `tests/`
directories, and `ScoreBGE` has pytest tests in [ScoreBGE/tests](ScoreBGE/tests).
