# Utility WDLs

This directory contains a collection of miscellaneous WDLs useful for some small tasks. Check below for documentation on each.

* [AnnotateVCF](#annotatevcf)
* [CollectBenchmarkSucceeded](#collectbenchmarksucceeded)
* [CombineTables](#combinetables)
* [CreateIGVSession](#createigvsession)
* [DetectPCANovelties](#detectpcanovelties)
* [Dipcall](#dipcall)
* [DownsampleAndCollectCoverage](#downsampleandcollectcoverage)
* [ExtractSampleFromVCF](#extractsamplefromvcf)
* [IndexCramOrBam](#indexcramorbam)
* [IntervalList2Bed](#intervallist2bed)
* [MatchFingerprints](#matchfingerprints)
* [MergeSingleSampleMinimacVcfs](#mergesinglesampleminimacvcfs)
* [PRSQC](#prsqc)
* [RNAMetrics](#rnametrics)


## AnnotateVCF

### Summary

This WDL takes a VCF and adds various annotations based on user input. Optionally, users can include a truth VCF to do a `vcfeval` benchmarking comparison, and use the resulting TP, FP, etc. labels in the output. The truth dataset will also be annotated using the input configuration. The final output will be table(s) with rows given by variants, and columns the resulting annotations. The families of possible annotations are:
* bed region membership: for each bed file provided, each variant gets a binary flag annotation with the corresponding label depending on if it lies in the region.
* reference / GC context: using bedtools, gather some statistics about GC content in a window around the variant, or occurrences of a custom sequence motif.
* general GATK annotations: using `VariantAnnotator`, add any of the annotations available in GATK, e.g. `Coverage`, or `FisherStrand`; requires reads be provided.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/AnnotateVCF.html) · [open locally](../../docs/viz/Utilities/WDLs/AnnotateVCF.html)

### Inputs
* `query_vcf`: VCF to annotate
* `query_vcf_index`: index for `query_vcf`
* `truth_vcf`: (optional) truth VCF to use for benchmarking
* `truth_vcf_index`: index for `truth_vcf`
* `truth_bed`: (optional) bed file for truth VCF evaluation regions
* `ref_fasta`: reference FASTA
* `ref_fasta_index`: index for `ref_fasta`
* `fasta_dict`: sequence dictionary for `ref_fasta`
* `bed_files`: list of bed files to use for region membership annotations
* `bed_labels`: list of labels for the bed files
* `add_gc_content`: (default: `true`) toggle to add GC content annotations
* `window_size`: (default: `25`) size of window to use for GC content annotations
* `include_N_count`: (default: `false`) toggle to include count of N's in GC content annotations
* `annotate_sequence`: (default: `false`) toggle to add reference sequence of `window_size` as an annotation
* `sequence_motif`: (optional) if provided, the number of occurrences of this sequence motif in the `window_size` will be added as an annotation
* `extra_query_fields`: (default: `[]`) list of extra fields to add to final table already in the query VCF
* `extra_truth_fields`: (default: `[]`) list of extra fields to add to final table already in the truth VCF
* `gatk_query_annotations`: (default: `[]`) list of GATK annotations to add to the query VCF
* `gatk_query_annotation_labels`: (default: `[]`) list of labels for the GATK annotations in the query VCF
* `gatk_truth_annotations`: (default: `[]`) list of GATK annotations to add to the truth VCF
* `gatk_truth_annotation_labels`: (default: `[]`) list of labels for the GATK annotations in the truth VCF
* `query_bam`: (optional) BAM file to use for GATK annotations
* `query_bam_index`: (optional) index for `query_bam`
* `gatk_jar`: (optional) path to GATK jar file if using custom build for user-defined `VariantAnnotation` classes; if not provided, the WDL will use the `gatk` command

The output tables will be split into SNPs and INDELs and will be in the format of a (gzipped) TSV file with a header row. This makes it easier for downstream analysis, like training separate models for SNPs and INDELs (see the `TrainEBMVariantAnnotations` notebook for some ideas).

## CollectBenchmarkSucceeded

### Summary

When running the `FindSamplesAndBenchmark` wdl with many samples, it sometimes happens that a few fail in Terra while most succeed.
Unfortunately, this means that the outputs of benchmarking the successful ones don't get compiled into one convenient .csv file to use for data analysis.
If you don't mind sacrificing the few that failed, or want to get started analyzing the successful ones ASAP, this wdl will automatically collect
the successful outputs and aggregate them into one .csv, similar to the last task of the benchmarking wdl.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/CollectBenchmarkSucceeded.html) · [open locally](../../docs/viz/Utilities/WDLs/CollectBenchmarkSucceeded.html)

### Inputs

* `namespace`: the first personalized part of your workspace URL; e.g. if you see `<my_project>/<my_workspace>` at the top
  in Terra, then this should be `<my_project>` as a string.
* `workspace_name`: specific name for your workspace, e.g. `<my_workspace>` in the last example.
* `submission_id`: the submission id for the `FindSamplesAndBenchmark` run, found from the "Job History" tab.


## CombineTables

### Summary

This WDL concatenates a list of tab-separated tables into one combined TSV using `pandas`. All input tables are read with
`#` treated as the comment character and then stacked row-wise, so the tables are expected to share a compatible set of
column headers. Optionally, you can tack on extra constant-valued columns to the combined output, which is useful for
labeling the rows of the result with metadata (e.g. a sample name or experiment tag) when aggregating scattered outputs
from another workflow. Runs on the `us.gcr.io/broad-dsde-methods/python-data-slim:1.0` docker.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/CombineTables.html) · [open locally](../../docs/viz/Utilities/WDLs/CombineTables.html)

### Inputs
* `tables`: list of TSV files to concatenate
* `output_name`: (default: `"combined_table"`) basename for the output file, which will be written as `<output_name>.tsv`
* `extra_column_names`: (default: `[]`) names of extra columns to add to the output
* `extra_column_values`: (default: `[]`) values to fill the extra columns with; one value per name in `extra_column_names`

### Outputs
* `combined_table`: the concatenated TSV, with any extra columns appended


## CreateIGVSession

### Summary

This workflow takes in optional lists of BAMs, VCFs, and interval files (`.interval_list` or `.bed`) and combines them together
into an IGV session .xml file. A reference must be provided, either by a hardcoded string ("hg38" or "hg19"), or by providing a 
path to the desired fasta. Input files are interpreted as WDL strings, so no localization occurs. Bucket paths are output in the .xml
session, so IGV will stream them directly from the cloud. This task is useful to add to the end of workflows that output lots of files
you might want to visualize together for analysis or debugging.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/CreateIGVSession.html) · [open locally](../../docs/viz/Utilities/WDLs/CreateIGVSession.html)

### Inputs

* `bams`: (optional) list of BAMs/CRAMs to add to session.
* `vcfs`: (optional) list of VCFs to add to session.
* `interval_files`: (optional) list of `.interval_list` or `.bed` files to add to session.
* `reference`: reference to use in IGV; must be either a `.fasta` file or one of the values: "hg38" or "hg19".
* `output_name`: (default = "igv_session") name for the output .xml file.


## DetectPCANovelties

### Summary

This WDL flags "novelties" (outliers) in a 2D PCA plot by building a concave hull around a training set and then checking
which test samples fall outside of it. The hull is an [alphashape](https://github.com/bellockk/alphashape), which is a
generalization of a convex hull whose tightness is controlled by the `alpha` parameter: larger values give a tighter
boundary around the training points. The workflow first calls `GenerateAlphashape` to fit the shape on the training data
and pickle it, then calls `DetectPCANoveltiesTask` to test each sample. A sample passes if it lies inside the shape, or
if its distance to the shape is below `distanceThreshold`; otherwise it is flagged as a novelty. A scatter plot showing
the shape, the training points, and the pass/fail labeled test points is also produced.

Both inputs are TSVs with a header row. The training file must have `PC1` and `PC2` columns, and the test file must have
`SAMPLE_ID`, `PC1`, and `PC2` columns. The underlying scripts live in this repo under
`Utilities/Dockers/Alphashape/` (`generate_alphashape.py` and `pca_novelty_detection.py`), and both tasks run on the
`us.gcr.io/broad-dsde-methods/kockan/alphashape` docker (pinned by digest in the WDL). This workflow is registered on
Dockstore as `DetectPCANovelties`. See also [PRSQC](#prsqc), which runs an inlined variant of this novelty check against
a pre-computed alphashape.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/DetectPCANovelties.html) · [open locally](../../docs/viz/Utilities/WDLs/DetectPCANovelties.html)

### Inputs
* `test`: TSV of samples to test, with columns `SAMPLE_ID`, `PC1`, `PC2`
* `training`: TSV of baseline/training samples, with columns `PC1`, `PC2`; used both to fit the alphashape and to draw the baseline points in the output plot
* `alpha`: (default: `8.0`) tightness of the alphashape fit; exposed on the `GenerateAlphashape` task
* `distanceThreshold`: (default: `0.01`) test points outside the shape but within this distance of it still pass; exposed on the `DetectPCANoveltiesTask` task

### Outputs
* `testSetPredictions`: a two-column TSV (no header) with one row per test sample, giving the sample id and either `PASS` or `FAIL`
* `runVisualization`: a PNG scatter plot of the alphashape with training points and pass/fail labeled test points


## Dipcall

### Summary

This WDL is a modified version of the Dockstore version of the [Dipcall](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/QC/wdl/tasks/dipcall.wdl)
pipeline. This workflow takes in a diploid assembly and calls variants, creating a VCF against your chosen reference.
This modified version allows for you to specify custom PAR regions for your reference so you can call haploid variants
when appropriate. Some data cleaning and indexing of the output VCF is also performed.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/Dipcall.html) · [open locally](../../docs/viz/Utilities/WDLs/Dipcall.html)

### Inputs

* `assemblyFastaPat`: the haploid assembly fasta for paternally inherited chromosomes
* `assemblyFastaMat`: the haploid assembly fasta for maternally inherited chromosomes
* `referenceFasta`: the reference to call variants against
* `isMaleSample`: set true if you would like to make haploid calls on X/Y outside the PAR region
* `custom_PAR_bed`: bed file denoting pseudoautosomal (PAR) regions for your reference
* `sample_name`: name to put for your sample in final output VCF
* `referenceIsHS38`: set true (default) if using hg38 reference


## DownsampleAndCollectCoverage

### Summary

The idea of this WDL is to do everything you need for a standard downsampling experiment. It takes in either CRAM or BAM files and downsamples them either according to a defined downsampling ratio or to a desired target coverage. If no downsampling ratio is defined then it will run `ColectWgsMetrics` to get the original mean coverage and determine the downsampling ratio based on that coverage and the desired target coverage. After downsampling using `DownsampleSam` the workflow will run `CollectWgsMetrics` once more and output the mean coverage of the downsampled CRAM (or BAM) file. This provides feedback with respect to the target coverage, because downsampling is always associated with some uncertainty. If `fail_if_below_coverage` is set, the workflow will fail if that downsampled mean coverage is below the provided threshold.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/DownsampleAndCollectCoverage.html) · [open locally](../../docs/viz/Utilities/WDLs/DownsampleAndCollectCoverage.html)

### Inputs
* `File input_cram`: Input BAM or CRAM
* `File input_cram_index`: Index file for input BAM or CRAM
* `File ref_fasta`: Reference FASTA
* `File ref_fasta_index`: Reference FASTA index
* `Float? downsample_probability`: Downsampling ratio. If not provided, the ratio will be determined based on the `target_coverage`.
* `Float? fail_if_below_coverage`: Fail the workflow if the downsampled mean coverage is below this value.
* `Float? target_coverage`: Target mean coverage for the downsampled CRAM file. **In order to use this input, do not provide `downsample_probability`, otherwise, that value will be used for downsampling.**
* `File? coverage_intervals`: If provided, the output downsampled mean coverage will be calculated based on these intervals. Additionally, these intervals will be used to calculate the original coverage if `target_coverage` is used.
* `String downsample_strategy = "ConstantMemory"`: See [DownsampleSam documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832708637467-DownsampleSam-Picard-).
* `Int read_length = 150`: See [CollectWgsMetrics documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832707851035-CollectWgsMetrics-Picard-).
* `Boolean use_fast_algorithm = true`: See [CollectWgsMetrics documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832707851035-CollectWgsMetrics-Picard-).
* `Boolean output_bam_instead_of_cram`: If set to true, the workflow will produce a downsampled output BAM (the default is CRAM).
* `String docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"`: Docker to use for both CollectWgsMetrics and DownsampleSam
* `File? picard_jar_override`: If provided, a Picard JAR file to use for both CollectWgsMetrics and DownsampleSam instead of the `gatk` command.
* `Int preemptible = 1`: Preemptible attempts

### Outputs
* `File downsampled_cram`: Downsampled CRAM
* `File downsampled_cram_index`: Downsampled CRAM index
* `Float downsampled_mean_coverage`: Mean coverage over the `coverage_intervals` (or the whole genome if not provided) for the downsampled CRAM file
* `File downsampled_wgs_metrics`: Output of CollectWgsMetrics run on the downsampled file
* `Float? original_mean_coverage`: The original mean coverage over the `coverage_intervals` (or the whole genome if not provided) of the input CRAM file, if `target_coverage` was used


## ExtractSampleFromVCF

### Summary

This WDL pulls a single sample out of a multi-sample callset VCF using GATK's
[SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/13832706016411-SelectVariants), writing a
block-gzipped VCF plus its index. Note the workflow itself is named `ExtractSingleSampleVCFFromCallset` (the file is
`ExtractSampleFromVCF.wdl`), which is also how it is registered on Dockstore. The task runs on the
`us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0` docker.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/ExtractSampleFromVCF.html) · [open locally](../../docs/viz/Utilities/WDLs/ExtractSampleFromVCF.html)

### Inputs
* `vcf`: the multi-sample callset VCF to extract from
* `basename`: basename for the output; the result is written as `<basename>.vcf.gz`
* `sampleName`: the sample to extract, passed to `SelectVariants --sample-name`

### Outputs
* `output_vcf`: the single-sample VCF
* `output_vcf_index`: the `.tbi` index for `output_vcf`


## IndexCramOrBam

### Summary 

Use this WDL to index a CRAM or BAM file, using `samtools`. The type is inferred using the file extension (either `.cram` or `.bam`). 

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/IndexCramOrBam.html) · [open locally](../../docs/viz/Utilities/WDLs/IndexCramOrBam.html)


## IntervalList2Bed

### Summary

This WDL takes in a list of interval files (either `.bed` or `.interval_list`) and converts the `.interval_list` files into `.bed`. The WDL checks if any of the provided files has a `.interval_list` extension, and then will call a conversion task on it if so. This means if all the files provided are `.bed`, then no tasks will be called, and the original list will be returned. This allows you to drop this task in to your workflows to extend pipeline functionality from accepting `.bed` inputs to also handle `.interval_list` files without penalizing users who provided `.bed` files with unnecessary extra tasks, which is ideal as many tools require specifically `.bed` lists.

If labels are provided, they will be returned in the new order of the `bed_files` output, which may be different than the originally given order. If labels are not provided, a list of `basename`s for the input files will be returned, in the correct order.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/IntervalList2Bed.html) · [open locally](../../docs/viz/Utilities/WDLs/IntervalList2Bed.html)

### Inputs
* `interval_files`: a list of `.bed` and/or `.interval_list` files
* `interval_labels`: an optional list of string labels to use for the corresponding interval file

### Outputs
* `bed_files`: a list of `.bed` files converted to the given inputs; note the order may have changed from the given list
* `bed_labels`: a list of labels for the `.bed` files corresponding to the user inputs, or the basename of the input files if the user did not provide any labels. Note the order may have changed, but the position of a label corresponds to the position of the file in `bed_files`.


## MatchFingerprints

### Summary

This WDL allows you to check fingerprints across two sets of files, and match them. There is an option to fail if files don't match fingerprints, allowing you to use this as a safety check on workflows that have paired files that must have matching samples. Alternatively, the WDL also has functionality to support finding matches across two batches, which can then be used downstream.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/MatchFingerprints.html) · [open locally](../../docs/viz/Utilities/WDLs/MatchFingerprints.html)

### Inputs
* `input_files`: a list of files to check the fingerprints of
* `input_indices`: list of corresponding index files for `input_files`
* `reference_files`: a list of files to use as a baseline when comparing fingerprints
* `reference_indices`: list of corresponding index files for `reference_files`
* `haplotype_map`: a haplotype map file used for Picard's `CrosscheckFingerprints` tool; see the docs [here](https://gatk.broadinstitute.org/hc/en-us/articles/13832766699291-CrosscheckFingerprints-Picard)
* `check_all_file_pairs`: (default: `true`) fingerprints pairs across *all* `input_file` and `reference_file` pairs when toggled `true`; otherwise fingerprints are only checked across files with the same index, and input lists must have the same length
* `fail_on_mismatch`: (default: `false`) toggle `true` to force your workflow to fail when fingerprinting fails to provide a "MATCH" for each comparison done; note for an individual comparison between `file1` and `file2`, if there are multiple samples/read groups/etc. being compared based on the mode selected, this check will pass (the workflow will NOT fail) if the resulting fingerprint summary file has *at least one* entry with a "MATCH"
* `check_only_matching_sample_names`: (default: `false`) toggle `true` to force the fingerprint comparison to have the same sample name across the files; required `crosscheck_by` be set to `SAMPLE`
* `crosscheck_by`: (default: `FILE`) controls at which level fingerprinting can happen; must be either `FILE`, `SAMPLE`, `LIBRARY`, or `READGROUP`
* `lod_threshold`: (default: -5) if the LOD fingerprinting score is less than this value, then the pair is a mismatch, and if it is greater than the negative of this value, then the pair is labeled a match

### Outputs
* `fingerprint_files`: a list of files output by `CrosscheckFingerprints` for each comparison made by the tool
* `matched_pairs`: a list of pairs of files that were detected to be matches using the set criteria; this list can be used/iterated over in other workflows to only act on pairs of files that are considered fingerprint matches


## MergeSingleSampleMinimacVcfs

### Summary

This WDL merges a large number of single-sample imputed VCFs produced by Minimac into a single multi-sample VCF. Because
every single-sample VCF from Minimac contains exactly the same sites in the same order, the merge can be done far more
cheaply than a general-purpose VCF merge: the `cut_paste_task` simply strips the headers, `paste`s the genotype columns
side by side, and re-attaches a single header. As a safety check, that task independently md5sums the site fields
(`CHROM`, `POS`, `ID`, `REF`, `ALT`, and `FORMAT`) of each input and fails if they do not all agree, so a mismatched set
of inputs cannot silently produce a corrupt merge.

The merge happens in two passes. First, the inputs are chunked into batches of `n_per_batch` and each batch is pasted
together. Then `interval_list` is split into `interval_scatter_count` pieces with Picard's `IntervalListTools`, each
batch VCF is subset to each interval piece with GATK `SelectVariants`, and the batches are pasted together per-interval
to produce the full sample set. Since `AF` and `R2` in the per-sample VCFs refer only to the imputation reference panel,
`reannotate_from_dosages` recomputes them across the merged cohort from the `DS` (dosage) field, replacing `AF` and `R2`
and dropping `INFO/MAF`. Finally the per-interval VCFs are concatenated with GATK `GatherVcfsCloud`.

In parallel, Hail's `sample_qc` is run on each batch to produce per-sample QC metrics, and those tables are concatenated
into one merged metrics TSV.

Note the Hail QC task initializes with `GRCh37` as the default reference. Dockers used: bcftools/bgzip
(`us.gcr.io/broad-dsde-methods/ckachulis/bcftools_bgzip`, pinned by digest) for the paste merge,
`us.gcr.io/broad-dsde-methods/samtools-suite:v1.1` for reannotation, `us.gcr.io/broad-dsde-methods/bcftools:v1.3` for
counting samples, `us.gcr.io/broad-gatk/gatk:4.3.0.0` for interval subsetting, `us.gcr.io/broad-gatk/gatk:4.5.0.0` for
the final gather, `us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039` for interval scattering,
`hailgenetics/hail:0.2.126-py3.11` for QC metrics, and `us.gcr.io/broad-dsde-methods/python-data-slim:1.1` for merging
the metrics. This workflow is registered on Dockstore and has an automated test configured in `test/watt_config.yml`.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/MergeSingleSampleMinimacVcfs.html) · [open locally](../../docs/viz/Utilities/WDLs/MergeSingleSampleMinimacVcfs.html)

### Inputs
* `vcfs`: list of single-sample imputed VCFs to merge; these must all contain the identical set of sites, in the same order
* `n_per_batch`: number of single-sample VCFs to paste together in each first-pass batch
* `output_basename`: basename for the merged VCF and QC metrics outputs
* `interval_list`: interval list used to scatter the second merge pass
* `interval_scatter_count`: (default: `100`) number of pieces to scatter `interval_list` into

### Outputs
* `merged_vcf`: the merged multi-sample VCF, with `AF` and `R2` recomputed from cohort dosages
* `merged_vcf_index`: the `.tbi` index for `merged_vcf`
* `merged_qc_metrics`: TSV of per-sample Hail `sample_qc` metrics across all samples


## PRSQC

### Summary

A simple QC workflow for polygenic risk scores, written for the PROGRESS VA project. It takes a table of PRS results and
runs two independent checks, then reports whether both passed.

`CheckScores` verifies that `prs_score`, `combined_risk_score`, `pc1`, and `pc2` all fall within the min/max bounds given
in the `acceptable_range` table, for every sample in the input. `DetectPCANovelties` checks that each sample's
(`pc1`, `pc2`) point falls inside a pre-computed [alphashape](https://github.com/bellockk/alphashape) (a concave hull
fit around a training population such as 1kG), or within `distance_threshold` of it, and emits a plot with each sample
colored green for pass or red for fail. Unlike the standalone [DetectPCANovelties](#detectpcanovelties) workflow, this
one takes the alphashape as an input rather than fitting it; you can generate one with
`Utilities/Dockers/Alphashape/generate_alphashape.py`. The default `alpha` used there is `8.0`, which was chosen
experimentally, and users training on something other than 1kG should be aware of that choice.

`CheckScores` runs on `us.gcr.io/broad-dsde-methods/python-data-slim:1.0` and `DetectPCANovelties` runs on the
`us.gcr.io/broad-dsde-methods/kockan/alphashape` docker (pinned by digest). This workflow is registered on Dockstore.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/PRSQC.html) · [open locally](../../docs/viz/Utilities/WDLs/PRSQC.html)

### Inputs
* `prs_full_risk`: TSV of PRS results indexed by `sample_id`, with columns `prs_score`, `combined_risk_score`, `pc1`, and `pc2`; works for either single-sample or multi-sample files
* `acceptable_range`: TSV indexed by metric name (`prs_score`, `combined_risk_score`, `pc1`, `pc2`) with `min` and `max` columns
* `output_basename`: basename for the output files
* `alphashape`: pickled alphashape bounding the (PC1, PC2) points of a training set
* `distance_threshold`: samples outside the alphashape but within this distance of it still pass

### Outputs
* `qc_passed`: `true` only if all metrics were within the acceptable range *and* all PCs were within the alphashape
* `pcs_within_shape`: `true` if every sample's (`pc1`, `pc2`) point passed the alphashape check
* `pca_qc_plot`: PNG of the alphashape with each sample plotted and colored by pass/fail


## RNAMetrics

### Summary

This WDL collects RNA-seq quality metrics for an aligned BAM by running two tools side by side on the same input:
[RNA-SeQC 2](https://github.com/getzlab/rnaseqc) (with a minimum mapping quality of 40) against a GTF annotation, and
Picard's [CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/13832708142235-CollectRnaSeqMetrics-Picard)
against a refFlat file and a ribosomal interval list. Note the Picard task hardcodes
`STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND`, so it assumes a dUTP-style stranded library. Output filenames are
derived from the input BAM's basename.

The two tasks live in `RNAMetricsTasks.wdl`, which `RNAMetrics.wdl` imports; if you run this outside of Dockstore make
sure that file sits alongside the main descriptor. Each task also runs a monitoring script from
`gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh` and returns its log. Dockers used:
`gcr.io/broad-cga-aarong-gtex/rnaseqc:latest` for RNA-SeQC 2 and `us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5` for
Picard. This workflow is registered on Dockstore as `RNAMetrics`.

**Interactive diagram:** [view on GitHub](https://raw.githack.com/broadinstitute/palantir-workflows/main/docs/viz/Utilities/WDLs/RNAMetrics.html) · [open locally](../../docs/viz/Utilities/WDLs/RNAMetrics.html)

### Inputs
* `inputBAM`: aligned RNA-seq BAM
* `inputBAMIndex`: index for `inputBAM`
* `referenceAnnotation`: GTF annotation passed to RNA-SeQC 2
* `referenceGenome`: reference FASTA
* `referenceGenomeIndex`: index for `referenceGenome`
* `refFlat`: refFlat-format gene annotations for `CollectRnaSeqMetrics`
* `ribosomalIntervals`: interval list of ribosomal regions for `CollectRnaSeqMetrics`

### Outputs
* `exonCV`: RNA-SeQC 2 per-exon coefficient of variation TSV
* `exonReads`: RNA-SeQC 2 per-exon read counts (`.gct`)
* `geneFragments`: RNA-SeQC 2 per-gene fragment counts (`.gct`)
* `geneReads`: RNA-SeQC 2 per-gene read counts (`.gct`)
* `geneTPM`: RNA-SeQC 2 per-gene TPM (`.gct`)
* `metrics`: RNA-SeQC 2 summary metrics TSV
* `rnaseqcMonitoringLog`: resource monitoring log for the RNA-SeQC 2 task
* `rnaMetrics`: Picard `CollectRnaSeqMetrics` output
* `collectRNASeqMetricsMonitoringLog`: resource monitoring log for the Picard task
