# gCNV (Germline Copy Number Variant) Workflows

This directory contains WDLs for calling germline copy number variants (CNVs) with the GATK
[GermlineCNVCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037593411-GermlineCNVCaller) tool chain,
plus the downstream filtering, reporting, merging, and QC steps that the Palantir team layers on top of it.

The core calling WDLs ([cnv_common_tasks.wdl](cnv_common_tasks.wdl),
[cnv_germline_cohort_workflow.wdl](cnv_germline_cohort_workflow.wdl), and
[cnv_germline_case_workflow.wdl](cnv_germline_case_workflow.wdl)) are derived from the germline CNV WDLs distributed
with GATK (`scripts/cnv_wdl/germline` in the [gatk repository](https://github.com/broadinstitute/gatk)) — the task
names, structure, GATK argument defaults, and header comments still match those upstream files closely. They were
brought into this repo already-modified in a single commit (`02882bd`, "BGE gCNV"), so the git history here does not
record which upstream version they were forked from. The remaining WDLs in this directory
([single_sample_cnv_germline_case_filter_workflow.wdl](single_sample_cnv_germline_case_filter_workflow.wdl),
[cnv_calling_and_merge_for_fabric.wdl](cnv_calling_and_merge_for_fabric.wdl),
[cohort_cnv_calling_and_merge_for_fabric.wdl](cohort_cnv_calling_and_merge_for_fabric.wdl), and
[CNVControlEventsQC.wdl](CNVControlEventsQC.wdl)) are local additions with no upstream counterpart.

Local changes relative to the upstream GATK germline WDLs that are visible in the source include:
* the cohort workflow emits extra "path list" outputs (plain text files listing the GCS paths of the sharded call
  tars, ploidy calls, and genotyped VCFs) via the `WritePathList` / `WritePathMatrix` tasks;
* the case workflow takes `normal_bams`/`normal_bais` as `Array[String]` (coerced to `File` inside the scatter) so
  that reads are streamed rather than localized;
* a single-sample variant that collapses `GermlineCNVCaller` (case mode) and `PostprocessGermlineCNVCalls` into one
  task, skips the interval scatter entirely, and adds panel-of-normals frequency annotation, `bcftools` soft
  filtering, and event-count QC;
* Fabric-oriented wrappers that reformat the CNV VCF, merge it with a short-variant VCF, and render an HTML report.

### Cohort mode vs. case mode

GATK gCNV is a two-stage method:

* **Cohort mode** (`GermlineCNVCaller --run-mode COHORT`, `DetermineGermlineContigPloidy` without `--model`) fits a
  denoising model and a contig-ploidy model from a cohort of samples, *and* emits calls for those same samples. This
  is how you build a panel of normals / model to be reused later.
  [CNVGermlineCohortWorkflow](#cnvgermlinecohortworkflow) implements this.
* **Case mode** (`GermlineCNVCaller --run-mode CASE`, `DetermineGermlineContigPloidy --model`) takes the model tars
  produced by a cohort run and calls new samples against them, which is far cheaper per sample.
  [CNVGermlineCaseWorkflow](#cnvgermlinecaseworkflow) (many samples, sharded model) and
  [SingleSampleGCNVAndFilterVCFs](#singlesamplegcnvandfiltervcfs) (one sample, one model shard) implement this.

Because cohort mode shards the model by interval, a case-mode run must use exactly the same interval sharding: the
`gcnv_model_tars` array from a cohort run has to line up index-for-index with the shards produced by
`ScatterIntervals` on the cohort's `filtered_intervals`.

This directory contains the following WDLs:
- [CNVGermlineCohortWorkflow](#cnvgermlinecohortworkflow) (`cnv_germline_cohort_workflow.wdl`): build gCNV
  denoising + ploidy models from a cohort of normals and genotype that cohort.
- [CNVGermlineCaseWorkflow](#cnvgermlinecaseworkflow) (`cnv_germline_case_workflow.wdl`): call CNVs for a batch of
  samples against an existing (sharded) cohort model.
- [SingleSampleGCNVAndFilterVCFs](#singlesamplegcnvandfiltervcfs)
  (`single_sample_cnv_germline_case_filter_workflow.wdl`): optimized single-sample, single-shard case-mode calling
  followed by panel-frequency annotation, filtering, and QC.
- [CNVCallingAndMergeForFabric](#cnvcallingandmergeforfabric) (`cnv_calling_and_merge_for_fabric.wdl`): runs
  `SingleSampleGCNVAndFilterVCFs`, merges the CNV calls with a short-variant VCF for delivery to Fabric, and renders
  an HTML CNV event report.
- [CohortCNVCallingAndMergeForFabric](#cohortcnvcallingandmergeforfabric)
  (`cohort_cnv_calling_and_merge_for_fabric.wdl`): the same end product, but starting from cohort-mode calling with
  no pre-existing model, using the cohort itself as its own panel of normals.
- [CNVControlEventsQC](#cnvcontroleventsqc) (`CNVControlEventsQC.wdl`): compares CNV calls on a control sample
  against a set of expected ("truth") events.
- [cnv_common_tasks.wdl](#cnv_common_taskswdl): shared task library imported by the calling workflows; not a
  runnable workflow.

### Dockstore and testing

Registered in [.dockstore.yml](../.dockstore.yml):
[CNVCallingAndMergeForFabric](#cnvcallingandmergeforfabric),
[CohortCNVCallingAndMergeForFabric](#cohortcnvcallingandmergeforfabric), and
[CNVControlEventsQC](#cnvcontroleventsqc). The other WDLs in this directory are not registered; they are either
imported by those three or run directly through Cromwell.

Correspondingly, [test/watt_config.yml](../test/watt_config.yml) defines a `base` test for exactly those three
workflows (`test/CNVCallingAndMergeForFabric`, `test/CohortCNVAndMergeForFabric`, and `test/CNVControlEventsQC`).
There is no automated test for `CNVGermlineCohortWorkflow`, `CNVGermlineCaseWorkflow`, or
`SingleSampleGCNVAndFilterVCFs` on their own, though the first and third are exercised indirectly through the two
Fabric workflows.

### Docker images

* **GATK**: every task that runs a GATK tool takes the image as a required `gatk_docker` workflow input, so there is
  no default baked into the WDL. The test inputs use `broadinstitute/gatk:4.4.0.0`. A `gatk4_jar_override` input is
  threaded through the GATK tasks (setting `GATK_LOCAL_JAR`) if you need to run a custom jar inside that image.
* `us.gcr.io/broad-dsde-methods/samtools-suite:v1.1` — panel annotation / `bcftools` filtering / event-count QC
  (`ExtractPoNFreqAnnotateFilterAndQC`).
* `us.gcr.io/broad-dsde-methods/r-gcnv-viz@sha256:c92a9a26...` (pinned by digest) — R/rmarkdown CNV event report
  (`GCNVVisualzation`).
* `us.gcr.io/broad-dsde-methods/python-h5py@sha256:7f8d5965...` (pinned by digest) — low-GC dropout metric
  (`LowGCDropoutQC`).
* `us.gcr.io/broad-dsde-methods/python-data-slim:1.1` — control-events QC (`CNVControlEventsQCTask`).
* `us.gcr.io/broad-dsde-methods/pysam:v1.1` — the (currently unused) `ReformatGCNVForFabric` task.
* `us.gcr.io/broad-dsde-methods/python-data-slim` (unpinned) — the (currently unused) `ExtractPoNFreq` task.
* `python:latest` — `WritePathList` / `WritePathMatrix` in the cohort workflow.

### A note on the GATK hyperparameters

`DetermineGermlineContigPloidy` and `GermlineCNVCaller` expose a large number of model / inference
hyperparameters. All of them are optional at the workflow level and, when left unset, the WDL passes the GATK
default explicitly on the command line. Rather than repeat every one below, they are grouped and the GATK defaults
that the WDLs hard-code are given; for the meaning of each, see the GATK tool documentation for
[DetermineGermlineContigPloidy](https://gatk.broadinstitute.org/hc/en-us/articles/360037593411) and
[GermlineCNVCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037593411-GermlineCNVCaller).

---

## CNVGermlineCohortWorkflow

Source: [cnv_germline_cohort_workflow.wdl](cnv_germline_cohort_workflow.wdl)

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/cnv_germline_cohort_workflow.html) · [open locally](../docs/viz/gCNV/cnv_germline_cohort_workflow.html)

### Summary

Builds a GATK gCNV denoising model and contig-ploidy model from a cohort of normal samples, and produces genotyped
CNV calls for those same samples. Supports both WGS and WES (for WES, supply the target interval list and set
`bin_length` to 0 so no binning is done; for WGS, supply intervals covering the chromosomes of interest and let the
default binning apply).

The steps are:
1. `PreprocessIntervals` — pad and bin the input intervals.
2. `AnnotateIntervals` — optional (on by default), adds GC content and, if tracks are supplied, mappability and
   segmental-duplication annotations; enables explicit GC correction in the denoising model.
3. `CollectCounts` (`CollectReadCounts`) — scattered over samples.
4. `FilterIntervals` — drops intervals by GC/mappability/segdup content and by read-count outlier criteria.
5. `DetermineGermlineContigPloidyCohortMode` — fits the ploidy model and calls per-sample contig ploidy.
6. `ScatterIntervals` (`IntervalListTools`) — split the filtered intervals into shards of roughly
   `num_intervals_per_scatter` intervals.
7. `GermlineCNVCallerCohortMode` — scattered over interval shards; fits the denoising model and emits per-sample,
   per-shard call tars.
8. `PostprocessGermlineCNVCalls` — scattered over samples; merges shards into genotyped interval and segment VCFs,
   and writes a per-sample QC status.
9. `CollectModelQualityMetrics`, `ScatterPloidyCallsBySample`, and the `WritePathList`/`WritePathMatrix` tasks that
   write text files of output paths (convenient for feeding large output arrays into downstream Terra data tables).

Note that `-L` intervals accept anything compatible with the GATK `-L` argument, and `blacklist_intervals` anything
compatible with `-XL` (useful for excluding centromeres, etc.).

### Inputs

**Required**
- **File intervals**: intervals over which to call; padded and binned by `PreprocessIntervals`.
- **Array[File]+ normal_bams** / **Array[File]+ normal_bais**: the cohort's BAMs/CRAMs and their indices. (Reads are
  declared `localization_optional` in `CollectCounts`, so they are streamed on backends that support it.)
- **String cohort_entity_id**: prefix used for all model/call output filenames.
- **File contig_ploidy_priors**: contig ploidy prior table for `DetermineGermlineContigPloidy`.
- **Int num_intervals_per_scatter**: approximate number of intervals per `GermlineCNVCaller` shard. See the comment
  on `ScatterIntervals` in [cnv_common_tasks.wdl](cnv_common_tasks.wdl) — `IntervalListTools` may produce shards
  larger than requested, so inspect the result.
- **File ref_fasta** / **File ref_fasta_fai** / **File ref_fasta_dict**: reference and its index/dictionary.
- **String gatk_docker**: GATK docker image.
- **Int ref_copy_number_autosomal_contigs**: reference copy number for autosomes (`PostprocessGermlineCNVCalls`).
- **Int maximum_number_events_per_sample** / **Int maximum_number_pass_events_per_sample**: QC thresholds; a sample
  whose genotyped segments VCF has at least this many non-ref (or non-ref PASS) events is flagged
  `EXCESSIVE_NUMBER_OF_EVENTS` / `EXCESSIVE_NUMBER_OF_PASS_EVENTS` instead of `PASS`. This is reported, not fatal.

**General optional**
- **File? blacklist_intervals**: intervals excluded at the `PreprocessIntervals` step (`-XL`).
- **Boolean? do_explicit_gc_correction**: defaults to `true`; when true `AnnotateIntervals` is run and its output is
  passed to `FilterIntervals` and `GermlineCNVCaller`.
- **File? gatk4_jar_override**: custom GATK jar (sets `GATK_LOCAL_JAR`, default `/root/gatk.jar`).
- **Int? preemptible_attempts**: applied to every task; each task has its own default (2 for the ploidy/caller
  tasks, 5 elsewhere).
- **String? gcs_project_for_requester_pays**: required if the BAM/CRAM lives in a requester-pays bucket.

**PreprocessIntervals**
- **Int? padding**: interval padding, default `250`.
- **Int? bin_length**: bin size, default `1000`; set to `0` to disable binning (typical for WES).

**AnnotateIntervals**
- **File? mappability_track_bed** / **File? mappability_track_bed_idx**
- **File? segmental_duplication_track_bed** / **File? segmental_duplication_track_bed_idx**
- **Int? feature_query_lookahead**: default `1000000`.
- **Int? mem_gb_for_annotate_intervals**: default 2 GB.

**FilterIntervals**
- **File? blacklist_intervals_for_filter_intervals**: `-XL` applied at the filtering step (separate from the
  `PreprocessIntervals` blacklist).
- **Float? minimum_gc_content** (`0.1`), **Float? maximum_gc_content** (`0.9`), **Float? minimum_mappability**
  (`0.9`), **Float? maximum_mappability** (`1.0`), **Float? minimum_segmental_duplication_content** (`0.0`),
  **Float? maximum_segmental_duplication_content** (`0.5`): annotation-based interval filters (only meaningful for
  annotations that were actually computed).
- **Int? low_count_filter_count_threshold** (`10`), **Float? low_count_filter_percentage_of_samples** (`50.0`),
  **Float? extreme_count_filter_minimum_percentile** (`1.0`), **Float? extreme_count_filter_maximum_percentile**
  (`99.0`), **Float? extreme_count_filter_percentage_of_samples** (`90.0`): count-based interval filters.
- **Int? mem_gb_for_filter_intervals**: default 7 GB.

**CollectCounts**
- **Array[String]? disabled_read_filters_for_collect_counts**: passed as `--disable-read-filter`.
- **String? collect_counts_format**: `HDF5` (default), `TSV`, or `TSV_GZ`. Anything else fails the task.
- **Boolean? collect_counts_enable_indexing**: default `false`; incompatible with `HDF5`.
- **Int? mem_gb_for_collect_counts**: default 7 GB.

**DetermineGermlineContigPloidy (cohort mode)**
- **Float? ploidy_mean_bias_standard_deviation** (`1`), **Float? ploidy_mapping_error_rate** (`0.3`),
  **Float? ploidy_global_psi_scale** (`0.001`), **Float? ploidy_sample_psi_scale** (`0.0001`).
- **Int? mem_gb_for_determine_germline_contig_ploidy** (default 7 GB) / **Int? cpu_for_determine_germline_contig_ploidy**
  (default 8; also sets `MKL_NUM_THREADS` and `OMP_NUM_THREADS`).

**GermlineCNVCaller (cohort mode)**

Caller parameters: **Float? gcnv_p_alt** (`5e-4`), **Float? gcnv_p_active** (`1e-1`),
**Float? gcnv_cnv_coherence_length** (`10000.0`), **Float? gcnv_class_coherence_length** (`10000.0`),
**Int? gcnv_max_copy_number** (`5`).

Denoising model parameters: **Int? gcnv_max_bias_factors** (`6`), **Float? gcnv_mapping_error_rate** (`0.01`),
**Float? gcnv_interval_psi_scale** (`0.01`), **Float? gcnv_sample_psi_scale** (`0.01`),
**Float? gcnv_depth_correction_tau** (`10000.0`), **Float? gcnv_log_mean_bias_standard_deviation** (`0.1`),
**Float? gcnv_init_ard_rel_unexplained_variance** (`0.1`), **Int? gcnv_num_gc_bins** (`20`),
**Float? gcnv_gc_curve_standard_deviation** (`1.0`), **String? gcnv_copy_number_posterior_expectation_mode**
(`HYBRID`), **Boolean? gcnv_enable_bias_factors** (`true`), **Int? gcnv_active_class_padding_hybrid_mode**
(`50000`).

Hybrid ADVI parameters: **Float? gcnv_learning_rate** (`0.05`), **Float? gcnv_adamax_beta_1** (`0.9`),
**Float? gcnv_adamax_beta_2** (`0.99`), **Int? gcnv_log_emission_samples_per_round** (`50`),
**Float? gcnv_log_emission_sampling_median_rel_error** (`0.005`), **Int? gcnv_log_emission_sampling_rounds** (`10`),
**Int? gcnv_max_advi_iter_first_epoch** (`5000`), **Int? gcnv_max_advi_iter_subsequent_epochs** (`100`),
**Int? gcnv_min_training_epochs** (`10`), **Int? gcnv_max_training_epochs** (`100`),
**Float? gcnv_initial_temperature** (`2.0`), **Int? gcnv_num_thermal_advi_iters** (`2500`),
**Int? gcnv_convergence_snr_averaging_window** (`500`), **Float? gcnv_convergence_snr_trigger_threshold** (`0.1`),
**Int? gcnv_convergence_snr_countdown_window** (`10`), **Int? gcnv_max_calling_iters** (`10`),
**Float? gcnv_caller_update_convergence_threshold** (`0.001`), **Float? gcnv_caller_internal_admixing_rate**
(`0.75`), **Float? gcnv_caller_external_admixing_rate** (`1.00`), **Boolean? gcnv_disable_annealing** (`false`).

Runtime: **Int? mem_gb_for_germline_cnv_caller** (default 7 GB), **Int? cpu_for_germline_cnv_caller** (default 8),
**Int? disk_for_germline_cnv_caller** (default 150 GB).

**PostprocessGermlineCNVCalls**
- **Array[String]? allosomal_contigs**: contigs to treat as allosomal (`--allosomal-contig`); their reference copy
  number comes from the ploidy calls rather than `ref_copy_number_autosomal_contigs`.
- **Int? mem_gb_for_postprocess_germline_cnv_calls** (default 7 GB) /
  **Int? disk_space_gb_for_postprocess_germline_cnv_calls** (default 40 GB).

### Outputs

- **File preprocessed_intervals**: padded/binned interval list.
- **Array[String] read_counts_entity_ids**: sample IDs, derived from the BAM/CRAM basenames.
- **Array[File] read_counts**: per-sample `CollectReadCounts` output.
- **File? annotated_intervals**: GC (and mappability/segdup) annotations, if `do_explicit_gc_correction`.
- **File filtered_intervals**: the interval list the model was trained on. Needed for case-mode runs.
- **File contig_ploidy_model_tar**: ploidy model; input to case-mode workflows.
- **File contig_ploidy_calls_tar**: cohort ploidy calls.
- **File contig_ploidy_calls_tar_path_list**: text file containing the GCS path of `contig_ploidy_calls_tar`.
- **Array[File] sample_contig_ploidy_calls_tars**: ploidy calls split into one tar per sample.
- **Array[File] gcnv_model_tars**: the denoising model, one tar per interval shard; input to case-mode workflows.
- **Array[Array[File]] gcnv_calls_tars**: raw call tars, indexed `[shard][sample]`.
- **File gcnv_calls_tars_path_list**: TSV of the above paths.
- **Array[File] gcnv_tracking_tars**: per-shard ADVI convergence tracking output.
- **Array[File] genotyped_intervals_vcfs** / **genotyped_intervals_vcf_indexes** and
  **Array[File] genotyped_segments_vcfs** / **genotyped_segments_vcf_indexes**: per-sample genotyped VCFs, plus
  **File genotyped_intervals_vcfs_path_list**, **genotyped_intervals_vcf_indexes_path_list**,
  **genotyped_segments_vcfs_path_list**, and **genotyped_segments_vcf_indexes_path_list** text files of their paths.
- **Array[File] denoised_copy_ratios**: per-sample denoised copy ratio TSVs.
- **Array[File] sample_qc_status_files** / **Array[String] sample_qc_status_strings**: `PASS`,
  `EXCESSIVE_NUMBER_OF_EVENTS`, or `EXCESSIVE_NUMBER_OF_PASS_EVENTS` per sample.
- **File model_qc_status_file** / **String model_qc_string**: `PASS` or `ALL_PRINCIPAL_COMPONENTS_USED` — the latter
  means no ARD component was pruned in some shard, i.e. `gcnv_max_bias_factors` may be too small.
- **Array[File] calling_configs**, **Array[File] denoising_configs**, **Array[File] gcnvkernel_version**,
  **Array[File] sharded_interval_lists**: per-shard metadata emitted by `GermlineCNVCaller`, required by
  `PostprocessGermlineCNVCalls`.

---

## CNVGermlineCaseWorkflow

Source: [cnv_germline_case_workflow.wdl](cnv_germline_case_workflow.wdl)

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/cnv_germline_case_workflow.html) · [open locally](../docs/viz/gCNV/cnv_germline_case_workflow.html)

### Summary

Calls CNVs for a batch of samples in case mode against models produced by
[CNVGermlineCohortWorkflow](#cnvgermlinecohortworkflow). It re-runs `PreprocessIntervals` and `CollectCounts`, then
`DetermineGermlineContigPloidy` in case mode against `contig_ploidy_model_tar`, scatters the cohort's
`filtered_intervals` the same way the cohort run did, runs `GermlineCNVCaller --run-mode CASE` once per shard
against the matching `gcnv_model_tars[shard]`, and finally runs `PostprocessGermlineCNVCalls` per sample.

Unlike the cohort workflow, `normal_bams`/`normal_bais` are declared as `Array[String]` and coerced to `File` inside
the scatter, so Cromwell does not localize them up front.

Case mode exposes a smaller set of hyperparameters than cohort mode — the denoising model is already fit, so the
model-structure parameters (bias factors, GC curve, interval psi scale, etc.) are not settable here.

### Inputs

**Required**
- **File intervals**: intervals for `PreprocessIntervals` (should match what the cohort run used).
- **File filtered_intervals**: the `filtered_intervals` output of the cohort run; scattered to define the shards.
- **Array[String]+ normal_bams** / **Array[String]+ normal_bais**: paths to the case BAMs/CRAMs and indices.
- **File contig_ploidy_model_tar**: `contig_ploidy_model_tar` from the cohort run.
- **Array[File]+ gcnv_model_tars**: `gcnv_model_tars` from the cohort run; must be in shard order.
- **Int num_intervals_per_scatter**: must reproduce the cohort run's sharding.
- **File ref_fasta** / **File ref_fasta_fai** / **File ref_fasta_dict**, **String gatk_docker**.
- **Int ref_copy_number_autosomal_contigs**.
- **Int maximum_number_events_per_sample** / **Int maximum_number_pass_events_per_sample**: QC thresholds, as above.

**General optional**
- **File? blacklist_intervals**, **File? gatk4_jar_override**, **Int? preemptible_attempts**,
  **String? gcs_project_for_requester_pays**: as in the cohort workflow (here the ploidy and caller tasks default to
  5 preemptible attempts).

**PreprocessIntervals**
- **Int? padding** (`250`), **Int? bin_length** (`1000`).

**CollectCounts**
- **Array[String]? disabled_read_filters_for_collect_counts**, **String? collect_counts_format** (`HDF5`),
  **Boolean? collect_counts_enable_indexing** (`false`), **Int? mem_gb_for_collect_counts** (7 GB).

**DetermineGermlineContigPloidy (case mode)**
- **Float? ploidy_mapping_error_rate** (`0.3`), **Float? ploidy_sample_psi_scale** (`0.0001`).
- **Int? mem_gb_for_determine_germline_contig_ploidy** (7 GB), **Int? cpu_for_determine_germline_contig_ploidy**
  (8), **Int? disk_for_determine_germline_contig_ploidy** (150 GB).

**GermlineCNVCaller (case mode)**

Caller parameters: **Float? gcnv_p_alt** (`5e-4`), **Float? gcnv_cnv_coherence_length** (`10000.0`),
**Int? gcnv_max_copy_number** (`5`).

Denoising parameters still settable in case mode: **Float? gcnv_mapping_error_rate** (`0.01`),
**Float? gcnv_sample_psi_scale** (`0.01`), **Float? gcnv_depth_correction_tau** (`10000.0`),
**String? gcnv_copy_number_posterior_expectation_mode** (`HYBRID`),
**Int? gcnv_active_class_padding_hybrid_mode** (`50000`).

Hybrid ADVI parameters (same names and defaults as the cohort workflow): **Float? gcnv_learning_rate** (`0.05`),
**Float? gcnv_adamax_beta_1** (`0.9`), **Float? gcnv_adamax_beta_2** (`0.99`),
**Int? gcnv_log_emission_samples_per_round** (`50`), **Float? gcnv_log_emission_sampling_median_rel_error**
(`0.005`), **Int? gcnv_log_emission_sampling_rounds** (`10`), **Int? gcnv_max_advi_iter_first_epoch** (`5000`),
**Int? gcnv_max_advi_iter_subsequent_epochs** (`100`), **Int? gcnv_min_training_epochs** (`10`),
**Int? gcnv_max_training_epochs** (`100`), **Float? gcnv_initial_temperature** (`2.0`),
**Int? gcnv_num_thermal_advi_iters** (`2500`), **Int? gcnv_convergence_snr_averaging_window** (`500`),
**Float? gcnv_convergence_snr_trigger_threshold** (`0.1`), **Int? gcnv_convergence_snr_countdown_window** (`10`),
**Int? gcnv_max_calling_iters** (`10`), **Float? gcnv_caller_update_convergence_threshold** (`0.001`),
**Float? gcnv_caller_internal_admixing_rate** (`0.75`), **Float? gcnv_caller_external_admixing_rate** (`1.00`),
**Boolean? gcnv_disable_annealing** (`false`).

Runtime: **Int? mem_gb_for_germline_cnv_caller** (7 GB), **Int? cpu_for_germline_cnv_caller** (8),
**Int? disk_for_germline_cnv_caller** (150 GB).

**PostprocessGermlineCNVCalls**
- **Array[String]? allosomal_contigs**, **Int? mem_gb_for_postprocess_germline_cnv_calls** (7 GB),
  **Int? disk_space_gb_for_postprocess_germline_cnv_calls** (40 GB).

### Outputs

- **File preprocessed_intervals**
- **Array[String] read_counts_entity_id** / **Array[File] read_counts**
- **Array[File] sample_contig_ploidy_calls_tars**: case-mode ploidy calls, split per sample.
- **Array[Array[File]] gcnv_calls_tars** (indexed `[shard][sample]`) and **Array[File] gcnv_tracking_tars**
- **Array[File] genotyped_intervals_vcfs** / **genotyped_intervals_vcf_indexes**
- **Array[File] genotyped_segments_vcfs** / **genotyped_segments_vcf_indexes**
- **Array[File] qc_status_files** / **Array[String] qc_status_strings**
- **Array[File] denoised_copy_ratios**
- **Array[File] sharded_interval_list**: per-shard `interval_list.tsv` from `GermlineCNVCaller`.

---

## SingleSampleGCNVAndFilterVCFs

Source: [single_sample_cnv_germline_case_filter_workflow.wdl](single_sample_cnv_germline_case_filter_workflow.wdl)

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/single_sample_cnv_germline_case_filter_workflow.html) · [open locally](../docs/viz/gCNV/single_sample_cnv_germline_case_filter_workflow.html)

### Summary

A single-sample, cost-optimized case-mode workflow with filtering and QC. Compared to
[CNVGermlineCaseWorkflow](#cnvgermlinecaseworkflow) it:

* takes an already-preprocessed interval list (no `PreprocessIntervals` call) and a **single** `gcnv_model_tar`,
  i.e. it assumes the cohort model was built as one shard, so there is no `ScatterIntervals` / scatter-gather;
* fuses `GermlineCNVCaller --run-mode CASE` and `PostprocessGermlineCNVCalls` into a single task
  (`GermlineCNVCallerCaseModeAndPostProcess`), avoiding an extra VM and the tar round trip;
* adds `ExtractPoNFreqAnnotateFilterAndQC`, which is where most of the added value is.

`ExtractPoNFreqAnnotateFilterAndQC` computes, for each called segment, how many samples in a panel of normals carry
an overlapping event of the same type. Events are expanded onto the exons (intervals) they cover; a panel sample
counts toward the event if more than `overlap_thresh` of the event's exonic length is covered by an event of the
same SVTYPE in that panel sample. The resulting `PANEL_COUNT` and `PANEL_FREQ` (count divided by the number of
panel samples) are annotated onto the VCF with `bcftools annotate`. The panel sample that matches the case sample's
own sample ID is excluded from the panel, so the workflow can safely be run with the case sample's own cohort as
the panel. Then each `filter_expressions[i]` is applied as a `bcftools filter` soft filter named
`filter_names[i]`, and finally the passing/total non-ref event counts are compared against the QC thresholds.

### Inputs

**Required**
- **File preprocessed_intervals**: preprocessed/binned interval list (also used as the "exon" definition for the
  panel-frequency calculation).
- **File normal_bam** / **File normal_bai**: the case sample's BAM/CRAM and index.
- **File contig_ploidy_model_tar**: ploidy model from a cohort run.
- **File gcnv_model_tar**: the (single) denoising model shard from a cohort run.
- **Array[File]+ pon_genotyped_segments_vcfs**: panel-of-normals genotyped segment VCFs used for the
  `PANEL_COUNT`/`PANEL_FREQ` annotation.
- **File ref_fasta** / **File ref_fasta_fai** / **File ref_fasta_dict**, **String gatk_docker**.
- **Int ref_copy_number_autosomal_contigs**.
- **Int maximum_number_events_per_sample** / **Int maximum_number_pass_events_per_sample**: QC thresholds. Note
  these are applied here as "more than" thresholds (`n > max` fails), and are evaluated *after* filtering, unlike
  the equivalent check inside `PostprocessGermlineCNVCalls`.

**General optional**
- **File? gatk4_jar_override**, **Int? preemptible_attempts**, **String? gcs_project_for_requester_pays**.

**CollectCounts**
- **Array[String]? disabled_read_filters_for_collect_counts**, **String? collect_counts_format** (`HDF5`),
  **Boolean? collect_counts_enable_indexing** (`false`), **Int? mem_gb_for_collect_counts** (7 GB).

**DetermineGermlineContigPloidy (case mode)**
- **Float? ploidy_mapping_error_rate** (`0.3`), **Float? ploidy_sample_psi_scale** (`0.0001`),
  **Int? mem_gb_for_determine_germline_contig_ploidy** (7 GB), **Int? cpu_for_determine_germline_contig_ploidy**
  (8), **Int? disk_for_determine_germline_contig_ploidy** (150 GB).

**GermlineCNVCaller (case mode) + postprocessing**
- Caller: **Float? gcnv_p_alt** (`5e-4`), **Float? gcnv_cnv_coherence_length** (`10000.0`),
  **Int? gcnv_max_copy_number** (`5`).
- Denoising: **Float? gcnv_mapping_error_rate** (`0.01`), **Float? gcnv_sample_psi_scale** (`0.01`),
  **Float? gcnv_depth_correction_tau** (`10000.0`), **String? gcnv_copy_number_posterior_expectation_mode**
  (`HYBRID`), **Int? gcnv_active_class_padding_hybrid_mode** (`50000`).
- Hybrid ADVI: **Float? gcnv_learning_rate** (`0.05`), **Float? gcnv_adamax_beta_1** (`0.9`),
  **Float? gcnv_adamax_beta_2** (`0.99`), **Int? gcnv_log_emission_samples_per_round** (`50`),
  **Float? gcnv_log_emission_sampling_median_rel_error** (`0.005`), **Int? gcnv_log_emission_sampling_rounds**
  (`10`), **Int? gcnv_max_advi_iter_first_epoch** (`5000`), **Int? gcnv_max_advi_iter_subsequent_epochs** (`100`),
  **Int? gcnv_min_training_epochs** (`10`), **Int? gcnv_max_training_epochs** (`100`),
  **Float? gcnv_initial_temperature** (`2.0`), **Int? gcnv_num_thermal_advi_iters** (`2500`),
  **Int? gcnv_convergence_snr_averaging_window** (`500`), **Float? gcnv_convergence_snr_trigger_threshold** (`0.1`),
  **Int? gcnv_convergence_snr_countdown_window** (`10`), **Int? gcnv_max_calling_iters** (`10`),
  **Float? gcnv_caller_update_convergence_threshold** (`0.001`), **Float? gcnv_caller_internal_admixing_rate**
  (`0.75`), **Float? gcnv_caller_external_admixing_rate** (`1.00`), **Boolean? gcnv_disable_annealing** (`false`).
- **Array[String]? allosomal_contigs**.
- Runtime: **Int? mem_gb_for_germline_cnv_caller** (7 GB), **Int? cpu_for_germline_cnv_caller** (8),
  **Int? disk_for_germline_cnv_caller** (150 GB).

**ExtractPoNFreqAnnotateFilterAndQC**
- **Float? overlap_thresh**: fraction of an event's exonic length that a panel event must cover for the panel sample
  to be counted. Task default `0.5`.
- **Array[String]? filter_expressions**: `bcftools filter -e` expressions. Task default:
  ```
  (GT=="alt" | GT=="mis") & ((FMT/CN>1 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400))
  (GT=="alt" | GT=="mis") & (INFO/PANEL_COUNT>1)
  ```
- **Array[String]? filter_names**: soft-filter labels, one per expression. Task default `['LowQual','PanelOverlap']`.
- **Int? mem_gb_for_extract_pon_freq** (16 GB) / **Int? disk_for_extract_pon_freq** (100 GB).

### Outputs

- **File filtered_vcf** / **File filtered_vcf_index** / **File filtered_vcf_md5sum**: the genotyped segments VCF
  annotated with `PANEL_FREQ`/`PANEL_COUNT` and soft-filtered, named `<sample>.filtered.genotyped-segments.vcf.gz`.
- **String read_counts_entity_id** / **File read_counts**
- **File contig_ploidy_calls_tar**
- **File gcnv_call_tar** / **File gcnv_tracking_tar**
- **File genotyped_intervals_vcf** / **File genotyped_intervals_vcf_index**
- **File genotyped_segments_vcf** / **File genotyped_segments_vcf_index**: the unfiltered, unannotated calls.
- **String qc_status_string**: `PASS`, `EXCESSIVE_NUMBER_OF_EVENTS`, or `EXCESSIVE_NUMBER_OF_PASS_EVENTS`.
- **Boolean qc_passed**: `qc_status_string == "PASS"`.
- **File cnv_metrics**: TSV with `sample`, `total_cnv_events`, `passing_cnv_events`.
- **File denoised_copy_ratios**
- **File interval_list**: the `interval_list.tsv` written by `GermlineCNVCaller`.

> The file also defines a standalone `ExtractPoNFreq` task that computes only the annotation TSV. It is not called
> by any workflow and appears to be superseded by `ExtractPoNFreqAnnotateFilterAndQC`.

---

## CNVCallingAndMergeForFabric

Source: [cnv_calling_and_merge_for_fabric.wdl](cnv_calling_and_merge_for_fabric.wdl) — registered in Dockstore,
tested by `test/CNVCallingAndMergeForFabric`.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/cnv_calling_and_merge_for_fabric.html) · [open locally](../docs/viz/gCNV/cnv_calling_and_merge_for_fabric.html)

### Summary

The single-sample delivery workflow. It runs
[SingleSampleGCNVAndFilterVCFs](#singlesamplegcnvandfiltervcfs), then:

* `ReformatAndMergeForFabric`: drops non-PASS records, rewrites `<DUP>` genotypes so the alt allele is on the last
  allele index and the remaining alleles are no-call (Fabric does not accept `./.` genotypes for these), adds a
  `CN` INFO header line, and merges the result with the sample's short-variant VCF using GATK `MergeVcfs` (with
  `-Dsamjdk.create_md5=true`). The task asserts that the CNV VCF basename and the short-variant VCF basename agree
  after stripping `.filtered.genotyped-segments.vcf.gz` and `.hard-filtered.vcf.gz` respectively, and fails if they
  do not.
* `GCNVVisualzation` (sic): renders an R Markdown HTML report with a genome-wide denoised copy-ratio panel, a
  copy-ratio-vs-GC-content panel, and, for each passing event, denoised copy ratio and panel-adjusted read count
  plots of the case sample against the panel over an automatically chosen window around the event.

Inputs not listed below (all the gCNV hyperparameters) are not exposed at this level, but can still be set through
Cromwell/Terra sub-workflow input syntax, e.g.
`CNVCallingAndMergeForFabric.SingleSampleGCNVAndFilterVCFs.filter_expressions` — the test inputs JSON does exactly
this.

### Inputs

- **File normal_bam** / **File normal_bai**: case sample reads.
- **File short_variant_vcf**: the sample's short-variant VCF to merge with. Must be named
  `<sample>.hard-filtered.vcf.gz` where `<sample>` matches the CNV VCF basename.
- **File contig_ploidy_model_tar**, **File preprocessed_intervals**, **File gcnv_model_tar**: the cohort model.
- **Array[File]+ gcnv_panel_genotyped_segments**: panel genotyped segment VCFs for the panel-frequency annotation.
- **Array[File]+ gcnv_panel_copy_ratios**: panel denoised copy ratio TSVs, for the report.
- **Array[File]+ gcnv_panel_read_counts**: panel `CollectReadCounts` HDF5s, for the report.
- **Float overlap_thresh**: default `0.5`; passed through to the panel-frequency calculation.
- **String gatk_docker**
- **Int maximum_number_events_per_sample** / **Int maximum_number_pass_events_per_sample**
- **Int ref_copy_number_autosomal_contigs**
- **File ref_fasta** / **File ref_fasta_fai** / **File ref_fasta_dict**
- **Array[String] allosomal_contigs**

### Outputs

- **File filtered_cnv_genotyped_segments_vcf** / **_index** / **_md5sum**: the annotated, soft-filtered CNV VCF.
- **File merged_vcf** / **File merged_vcf_index** / **File merged_vcf_md5sum**: short variants + PASS CNVs, named
  `<sample>.merged.vcf.gz`. Produced regardless of QC status.
- **Boolean qc_passed**
- **File cnv_metrics**: the per-sample event-count TSV.
- **File cnv_event_report**: the HTML CNV event report.

> The file also defines `ReformatGCNVForFabric` and `MergeVcfs` tasks, which together do what
> `ReformatAndMergeForFabric` now does in one task. Neither is called by any workflow.

---

## CohortCNVCallingAndMergeForFabric

Source: [cohort_cnv_calling_and_merge_for_fabric.wdl](cohort_cnv_calling_and_merge_for_fabric.wdl) — registered in
Dockstore, tested by `test/CohortCNVAndMergeForFabric`.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/cohort_cnv_calling_and_merge_for_fabric.html) · [open locally](../docs/viz/gCNV/cohort_cnv_calling_and_merge_for_fabric.html)

### Summary

Produces the same per-sample deliverables as [CNVCallingAndMergeForFabric](#cnvcallingandmergeforfabric), but for a
cohort with no pre-existing gCNV model. It runs [CNVGermlineCohortWorkflow](#cnvgermlinecohortworkflow) to build the
model and genotype the cohort, then scatters over the cohort's samples and, for each one, reuses tasks from the
single-sample path:

* `ExtractPoNFreqAnnotateFilterAndQC` with the whole cohort as the panel (the sample's own VCF is dropped from the
  panel inside the task by sample ID);
* `ReformatAndMergeForFabric` against the corresponding entry of `short_variant_vcfs` (**the arrays are matched by
  index**, so `short_variant_vcfs` must be in the same order as `normal_bams`);
* `GCNVVisualzation` with the cohort as the panel;
* `LowGCDropoutQC`, a cohort-workflow-only metric: it reads the sample's `CollectReadCounts` HDF5 and the
  `annotated_intervals` GC content, restricts to bins with GC content in (0.25, 0.30), and reports the fraction of
  those bins whose count is below 25% of the mean count in that GC band. This requires the cohort run's
  `annotated_intervals`, i.e. `do_explicit_gc_correction` must be left at its default of true.

Only a subset of the cohort workflow's inputs is exposed at this level; the rest can be reached with
`CohortCNVCallingAndMergeForFabric.CNVGermlineCohortWorkflow.<input>` sub-workflow syntax (the test inputs JSON does
this for the memory settings).

### Inputs

**Required**
- **File intervals**: intervals to call over.
- **Array[File]+ normal_bams** / **Array[File]+ normal_bais**
- **String cohort_entity_id**
- **File contig_ploidy_priors**
- **Int num_intervals_per_scatter**
- **File ref_fasta** / **File ref_fasta_fai** / **File ref_fasta_dict**
- **String gatk_docker**
- **Array[File] short_variant_vcfs**: one per sample, index-matched to `normal_bams`.
- **Int ref_copy_number_autosomal_contigs**
- **Int maximum_number_events_per_sample** / **Int maximum_number_pass_events_per_sample**: used both for the
  cohort workflow's internal QC and for the post-filtering QC.

**Optional**
- **Int? padding** / **Int? bin_length**: `PreprocessIntervals` settings.
- **Int? mem_gb_for_germline_cnv_caller** / **Int? cpu_for_germline_cnv_caller**
- **Int? mem_gb_for_postprocess_germline_cnv_calls** / **Int? disk_space_gb_for_postprocess_germline_cnv_calls**
- **Array[String] filter_expressions**: default
  ```
  (GT=="alt" | GT=="mis") & ((FMT/CN>1 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400))
  (GT=="alt" | GT=="mis") & (INFO/PANEL_FREQ>0.02)
  ```
  Note this differs from the single-sample workflow's default, which uses `INFO/PANEL_COUNT>1` — a frequency
  threshold makes more sense when the panel size is the cohort size.
- **Array[String] filter_names**: default `['LowQual','PanelOverlap']`.
- **Int? mem_gb_for_extract_pon_freq** / **Int? disk_for_extract_pon_freq**
- **Float? overlap_thresh**: panel overlap fraction, task default `0.5`.

### Outputs

All outputs are arrays, one entry per sample, in cohort order.

- **Array[File] filtered_cnv_genotyped_segments_vcf** / **_index** / **_md5sum**
- **Array[File] merged_vcf** / **merged_vcf_index** / **merged_vcf_md5sum**
- **Array[Boolean] qc_passed**
- **Array[File] cnv_metrics**
- **Array[File] cnv_event_report**
- **Array[File] low_gc_dropout_metric**: TSV with `sample` and `low_gc_dropout_frac`.

---

## CNVControlEventsQC

Source: [CNVControlEventsQC.wdl](CNVControlEventsQC.wdl) — registered in Dockstore, tested by
`test/CNVControlEventsQC`.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/gCNV/CNVControlEventsQC.html) · [open locally](../docs/viz/gCNV/CNVControlEventsQC.html)

### Summary

Process-control QC for a control sample (e.g. NA12878 or NA24385) run through the CNV pipeline. It compares the
PASS, non-ref events in a control sample's genotyped segments VCF against a TSV of expected events, using the same
exon-expansion logic as the panel-frequency calculation: both truth and eval events are expanded onto the exons they
overlap, and overlap is measured in exonic base pairs by contig, exon index, and SVTYPE.

If the eval VCF contains no passing events, the workflow does not fail; it reports zeros (and
`expected_events_not_seen` equal to the number of truth events).

The truth TSV is produced by [CNVControlEventsQCGenerateTruth.ipynb](CNVControlEventsQCGenerateTruth.ipynb), a
notebook that takes the filtered genotyped-segments VCFs from several pipeline runs of the same control sample,
expands them onto exons, keeps exon/SVTYPE combinations seen in at least 4 of the runs, and merges runs of
consecutive exons back into events. The notebook contains the (hard-coded) GCS paths of the runs used to build the
current NA12878 and NA24385 truth sets. The `exon_intervals` input to this workflow should be the same interval
list used to generate the truth.

### Inputs

- **File eval_control_sample**: genotyped segments VCF for the control sample to evaluate. It is read with
  `pandas.read_csv` without an explicit `compression` argument, so compression is inferred from the file extension;
  the test input is an uncompressed `.vcf`.
- **File control_sample_common_events**: the truth TSV (columns `contig`, `svtype`, `start`, `end`) generated by the
  notebook.
- **File exon_intervals**: Picard-style interval list defining the exons/bins.

### Outputs

- **Float sensitivity**: exonic overlap bp between truth and eval, divided by total truth exonic bp.
- **Float precision**: exonic overlap bp between truth and eval, divided by total eval exonic bp.
- **Float expected_events_seen**: sum over truth events of the fraction of each event's exonic length covered by
  eval events — i.e. a fractional count of recovered truth events.
- **Float expected_events_not_seen**: number of truth events minus `expected_events_seen`.
- **Float unexpected_events_seen**: number of eval events minus the equivalent fractional count of eval events
  covered by truth.

---

## cnv_common_tasks.wdl

Source: [cnv_common_tasks.wdl](cnv_common_tasks.wdl)

This is a task library, not a runnable workflow. It is imported (as `CNVTasks`) by
[cnv_germline_cohort_workflow.wdl](cnv_germline_cohort_workflow.wdl),
[cnv_germline_case_workflow.wdl](cnv_germline_case_workflow.wdl), and
[single_sample_cnv_germline_case_filter_workflow.wdl](single_sample_cnv_germline_case_filter_workflow.wdl), and
holds the tasks that cohort mode and case mode share. Every task in it takes `gatk_docker` and the usual optional
`mem_gb` / `disk_space_gb` / `use_ssd` / `cpu` / `preemptible_attempts` runtime overrides plus an optional
`gatk4_jar_override`; the calling workflows expose only some of these, and the per-task defaults are noted in the
workflow sections above.

Tasks provided:

- **PreprocessIntervals**: GATK `PreprocessIntervals`. Pads (`padding`, default 250) and bins (`bin_length`,
  default 1000) the input intervals, optionally excluding `blacklist_intervals`. Output basename is derived from the
  input interval list (or `wgs` if no intervals given).
- **AnnotateIntervals**: GATK `AnnotateIntervals`. Adds GC content, and mappability / segmental-duplication
  annotations if the corresponding BED tracks are provided.
- **FilterIntervals**: GATK `FilterIntervals`. Filters bins on the annotations above and on read-count statistics
  across the cohort.
- **CollectCounts**: GATK `CollectReadCounts` (plus `bgzip` and `IndexFeatureFile` when requested). Supports
  `HDF5` (default), `TSV`, and `TSV_GZ` output; validates that the format is known and that indexing is not
  requested with HDF5. The sample/entity ID is derived from the BAM/CRAM basename, and the reads are declared
  `localization_optional`.
- **CollectAllelicCounts**: GATK `CollectAllelicCounts` over a `common_sites` list (`minimum_base_quality` default
  20). Not currently called by any workflow in this directory.
- **ScatterIntervals**: Picard/GATK `IntervalListTools` with `SUBDIVISION_MODE=INTERVAL_COUNT` and
  `SCATTER_CONTENT=num_intervals_per_scatter`, renaming outputs to sort correctly. Short-circuits to a plain copy
  when only one shard is needed. The header comment warns that `IntervalListTools` can produce shards larger than
  requested.
- **PostprocessGermlineCNVCalls**: GATK `PostprocessGermlineCNVCalls`. Unpacks the per-shard call and model tars for
  one sample, produces genotyped interval and segment VCFs and a denoised copy ratio TSV, and then writes a QC
  status (`PASS` / `EXCESSIVE_NUMBER_OF_EVENTS` / `EXCESSIVE_NUMBER_OF_PASS_EVENTS`) based on the non-ref event
  counts in the segments VCF. Also accepts optional `intervals_vcf` / `clustered_vcf` / reference inputs for the
  breakpoint-clustering arguments, which none of the workflows here currently set.
- **CollectModelQualityMetrics**: checks each shard's `mu_ard_u_log__.tsv` for at least one positive ARD component;
  emits `PASS` or `ALL_PRINCIPAL_COMPONENTS_USED`.
- **ScatterPloidyCallsBySample**: splits a cohort/case contig-ploidy calls tar into one tar per sample, named so
  that the glob returns them in sample order.
- **SplitInputArray**: reshapes a flat `Array[String]` into a 2-D array with `num_inputs_in_scatter_block` columns.
  Not currently called by any workflow in this directory.
