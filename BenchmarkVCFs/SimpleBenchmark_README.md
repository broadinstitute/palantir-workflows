# SimpleBenchmark Guide

This guide describes an overview of the [SimpleBenchmark](SimpleBenchmark.wdl) WDL used for evaluating VCF performance against a known baseline. This includes how to toggle inputs for different goals and a description of the code structure.

Some useful resources referenced below:
* `bcftools` [Filtering Expressions](https://samtools.github.io/bcftools/bcftools.html#expressions)
* [RTG manual](https://cdn.jsdelivr.net/gh/RealTimeGenomics/rtg-core@master/installer/resources/core/RTGOperationsManual.pdf),
particularly the section on `vcfeval`

---

## Acknowledgements

This tool was heavily inspired by the previously developed BenchmarkVCFs and FindSamplesAndBenchmark tools, with some code lifted directly or slightly modified from both. The overall layout of this workflow follows the GA4GH Benchmarking guidelines as outlined [here](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/README.md). 

---

## Table of Contents

1. [Usage](#usage)
   1. [Workflow Inputs](#workflow-inputs)
   2. [Interpreting Output Data](#interpreting-output-data)
2. [Structure of the Workflow](#structure-of-the-workflow)
   1. [Classify Variants](#classify-variants)
   2. [Compute Statistics](#compute-statistics)
   3. [Combine Statistics](#combine-statistics)
3. [Ploidy Defaults](#ploidy-defaults)


---

## Usage

### Workflow Inputs

The following are the top-level inputs for the workflow. We borrow the `vcfeval` terminology in referring to the paradigm
of a "base"line VCF and a "call"set VCF benchmarked against it.

* `base_vcf` - the VCF to use as baseline.
* `base_vcf_index` - the index of the VCF to use as baseline.
* `base_output_sample_name` - the name to appear in output data to label the baseline sample.
* `base_vcf_sample_name` - (optional) the name appearing in the VCF to label this sample; to be used only for extracting that sample from a multisample callset. **Tip**: If you expect to use the same baseline repeatedly, consider extracting it from a multisample callset once, and then use that as input here to save time repeatedly subsetting a large VCF.

* `query_vcf` - the VCF to be evaluated against the baseline.
* `query_vcf_index` - the index of the VCF to be evaluated.
* `query_output_sample_name` - the name to appear in output data to label the evaluated sample.
* `query_vcf_sample_name` - (optional) the name appearing in the VCF to label this sample; to be used only for extracting that sample from a multisample callset.

* `ref_fasta` - reference fasta file matching the input VCFs.
* `ref_index` - reference fasta index.

* `stratifier_intervals` - (default = `[]`) list of interval files to stratify by to compute benchmarking statistics over. These can be any of the allowed interval file formats supported by GATK, including bed and Picard interval_list files.**Regardless of input, an "empty" interval file will be added to perform benchmarking over without any subseting, which will have a corresponding empty label.** This is important when computing summary statistics in groups, since some programs(like pandas groupby) default to excluding entries with null values. The naming convention for the empty file is left empty to avoid collisions with user input names on other interval files.
* `stratifier_labels` - (default = `[]`) list of names for the given `strat_intervals`. The order must match, and the sizes must agree.

* `evaluation_intervals` - if provided, only variants inside this bed file will be counted in statistics. Note this is slightly different from using a `strat_interval` because `vcfeval` will process variants on the boundary slightly differently here (see the RTG manual for details). This is intended to be used when one has a "high confidence" region for known truth samples.
* `score_field` - (default = `"GQ"`) an optional field to use for stratifying ROC data. The `vcfeval` will be run to produce the proper ROC data with the given string input to `--vcf-score-field`. Example values include `"GQ"`, `INFO.<NAME>`, and `FORMAT.<NAME>`. See the `vcfeval` documentation for more information **WARNING**: Some values might not make sense to stratify when interested in *both* precision and recall statistics. For example, with the value `FORMAT.DP` it might make sense to evaluate precision of calls made by depth in the callset, but often the baseline/truth doesn't have any intrinsic DP associated with each variant, and so any present in the baseline may lead to misleading recall statistics. See the section on [Interpreting Output Data](#interpreting-output-data) below for details.

* `experiment` - an optional string to populate "Experiment" column in final output tables. Useful for incorporating Terra table information about experimental groups to do grouped analyses downstream, or for easier integration into provided visualizer tools.
* `extra_column_name` - an optional string title for an additional column to add to final output tables, for downstream analysis.
* `extra_column_value` - values to populate in the optional `extra_column`, if specified. To make downstream analysis easier by incorporating Terra table information, e.g. `MEAN_COVERAGE`.

* `create_igv_session` - if true, an XML for an IGV session with relevant files will be created and output
* `optional_igv_bams` - an optional list of bams to include in the IGV session output; these are considered `String`s, so no overhead runtime cost is accrued from larger files.
* `igv_session_name` - the name of the XML file for IGV.

* `preemptible` - number of retries allowed on preemptible machines.

* `par_bed` - if provided (for VCFEval task), variants inside the X and Y choromosmes, which are also outside this bed file, will have the --squash-ploidy flag applied to them. Other regions will be treated as normal. See the default ploidy matcings [here](#ploidy-defaults).

In addition, there are some `vcfeval` options that can be toggled using the Boolean inputs `require_matching_genotypes`,`enable_ref_overlap`, and `passing_only`. Check the RTG manual for more details on how these can affect output. For most users, these can be safely left as default.

Each task also has inputs to control the amount of disk, CPU, and memory allocated. These can be toggled to optimize
performance on different inputs.


### Interpreting Output Data

There are four files output by the workflow. They are:

* `SimpleSummary` - This file contains counts of the usual TP, FP, etc. statistics, along with computed Precision, Recall, and F1_Score stats. These are broken down by user-provided intervals, and by zygosity of variants. 
* `IndelDistributionStats` - This file contains distributions of the TP, FP, etc. statistics by length of INDEL, where negative values correspond to deletions. For convenience, the column `INDEL_Type` labels the variant class as `Ins` or `Del` accordingly, so users can readily group and process their own statistics per class.
* `SNPSubstitutionStats` - This file contains distributions of the TP, FP, etc. statistics by SNP substitution
(e.g. `A>C`, `A>G`, etc.). For convenience, the column `Substitution_Type` labels the variant class as `Ti` or `Tv` depending on if the substitution is a transition or transversion accordingly, so users can readily group and process their own statistics per class.
* `ROCStats` - This file is the ROC outputs by `vcfeval`, relative to the whole genome. A `score_field` column records the user input score label for the given row.

We follow the same convention as RTG when computing Precision / Recall statistics, namely we use:

$$ Prec = \frac{TP_{call}}{(TP_{call} + FP)}, $$

$$ Recall = \frac{TP_{base}}{(TP_{base} + FN)}, $$

and

$$ F_1 = 2 \cdot \frac{Prec \cdot Recall}{Prec + Recall}.$$

---

## Structure of the Workflow

The workflow has the following main high-level tasks:

1. **Classify Variants** in both VCFs using `vcfeval`.
2. **Compute Statistics** on labeled variants in different subsets of types of variants in different regions.
3. **Combine Statistics** into a few final files for output.

Some details are provided for each step below.

### Classify Variants

We run `vcfeval` just once to perform haplotype-aware variant comparison. This is a sophisticated way to comparing variants across VCFs which accounts for potentially different methods of representing haplotype-equivalent variants. Each variant in the files can be categorized as true positive (TP), false positive (FP), and false negative (FN) based on how the call variants compare with the base. The command outputs a "combined" VCF file containing "consensus" calls between the two VCFs and classifier labels.

In addition, the command outputs ROC-style statistics split by user-input score fields (e.g. stratify TP, FP, FN counts by GQ scores, etc.). Some data cleaning is performed in this step on the ROC files output by `vcfeval`.

### Compute Statistics

In this step, we subset the combined VCF based on the TPs, FPs, and FNs and run `bcftools stats` on them to collect an array of statistics, particularly the count of SNPs and INDELs in each of these categories. This is performed once for each user-input region and across all, heterozygous, or homozygous variants, which must be communicated using the `bcftools` expressive filtering language. A large amount of data cleaning is also done in this step on the outputs to prepare for them to be compiled into the final outputs later.

### Combine Statistics

This step collects all the statistics files generated previously over the different interval files and combines them into one large file. Each row will contain unique identifiers to explain the interval file subset to and any variant selector conditions used, along with user-provided input sample names, so that users can safely concatenate outputs across different runs of the workflow for multi-sample evaluation analysis. There is also a little bit more cleaning to create a "simple summary" output, which is reported alongside the more in-depth statistics collected before.

---

## Ploidy Defaults

Different VCFs may choose to represent a given genotype differently. VCFeval considers some of these discrepancies to be matches, while others it considers to be mismatches. This is relevant for the X and Y chromosomes.

### Chart of Ploidy Deafults

This table shows the behavior of how `vcfeval` labels equal variants with differing `GT` fields. To reproduce the results in this table, see the `VCFEval_Ploidy` analysis done in the private hydro.gen repo.


| First  |Second| Result   |
| ------ | ---- |----------|
| 1/0    | 1    | Mismatch |
| 1/0    | 1/1  | Mismatch |
| 1/0    | 1/.  | Mismatch |
| 1/.    | 1    | Mismatch |
| 1/.    | 1/1  | Mismatch |
| 1/1    | 1    | Match    |
