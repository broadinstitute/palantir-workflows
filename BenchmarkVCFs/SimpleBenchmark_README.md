# SimpleBenchmark Guide

This guide describes an overview of the [SimpleBenchmark](SimpleBenchmark.wdl) WDL used for evaluating VCF performance 
against a known baseline. This includes how to toggle inputs for different goals, a description of the code structure, 
and comparison with previously developed tools.

Some useful resources referenced below:
* `bcftools` [Filtering Expressions](https://samtools.github.io/bcftools/bcftools.html#expressions)
* [RTG manual](https://cdn.jsdelivr.net/gh/RealTimeGenomics/rtg-core@master/installer/resources/core/RTGOperationsManual.pdf), 
particularly the section on `vcfeval`

---

## Acknowledgements

This tool was heavily inspired by the [BenchmarkVCFs](https://github.com/broadinstitute/palantir-workflows/blob/main/BenchmarkVCFs/BenchmarkVCFs.wdl) 
and [FindSamplesAndBenchmark](https://github.com/broadinstitute/palantir-workflows/blob/main/BenchmarkVCFs/FindSamplesAndBenchmark.wdl) tools, 
with some code lifted directly or slightly modified from both. Although this tool is mainly meant to supersede BenchmarkVCFs in most use cases, 
there are a small number of changes to the outputs which are summarized at the end of this article that old users will have to 
adapt to before sticking it into their pipelines.

---

## Table of Contents

1. [Usage](#usage)
    1. [Inputs](#inputs)
    2. [Interpreting Output Data](#interpreting-output-data)
        1. [Some Warnings on Interpreting Output Statistics](#some-warnings-on-interpreting-output-statistics)
2. [Structure of the Workflow](#structure-of-the-workflow)
    1. [Subset](#subset)
    2. [Classify Variants](#classify-variants)
    3. [Compute Statistics](#compute-statistics)
    4. [Combine Statistics](#combine-statistics)
3. [Comparison With BenchmarkVCFs](#comparison-with-benchmarkvcfs)

---

## Usage

### Inputs

The following are the top-level inputs for the workflow. We borrow the `vcfeval` terminology in referring to the paradigm 
of a "base"line VCF and a "call"set VCF benchmarked against it.

* `base_vcf` - the VCF to use as baseline.
* `base_vcf_index` - the index of the VCF to use as baseline.
* `base_output_sample_name` - the name to appear in output data to label the baseline sample.
* `base_vcf_sample_name` - (optional) the name appearing in the VCF to label this sample; to be used only for extracting that sample 
from a multisample callset. **Tip**: If you expect to use the same baseline repeatedly, consider extracting it from a 
multisample callset once, and then use that as input here to save time repeatedly subsetting a large VCF.

* `call_vcf` - the VCF ("callset") to be evaluated.
* `call_vcf_index` - the index of the VCF to be evaluated.
* `call_output_sample_name` - the name to appear in output data to label the evaluated sample.
* `call_vcf_sample_name` - (optional) the name appearing in the VCF to label this sample; to be used only for extracting that sample 
from a multisample callset.

* `ref_fasta` - reference fasta file matching the input VCFs.
* `ref_index` - reference fasta index.
* `ref_dict` - reference dictionary.

* `strat_intervals` - (default = `[]`) list of interval files to stratify by to compute benchmarking statistics over. 
These can be any of the allowed interval file formats supported by GATK, including bed and Picard interval_list files. 
**Regardless of input, an "empty" interval file will be added to perform benchmarking over without any subseting, which 
will have a corresponding empty label.** This is important when computing summary statistics in groups, since some programs 
(like pandas groupby) default to excluding entries with null values. The naming convention for the empty file is left empty 
to avoid collisions with user input names on other interval files.
* `strat_labels` - (default = `[]`) list of names for the given `strat_intervals`. The order must match, and the sizes must agree.
* `interval_padding` - (default = 0) an integer number of bases to pad the given `strat_intervals` by when subsetting. 
Useful for toggling to see how much statistics change when adding more bases, to help control for variants occuring near 
the boundaries of given intervals.
* `subset_gatk_tag` - (default = "4.2.3.0") a (string) version to use for GATK when performing subsetting.

* `bcf_selectors` - (default = `[]`) list of strings representing bcftools filtering expressions (see reference above). 
For each provided, `bcftools stats` will be run using the filter expression and counting TP, FP, etc. on the resulting 
subset of variants. Note that the selectors `["", "GT='het'", "GT='hom'"]` will be appended to your inputs to cover the 
default cases of no selection, Het variants, and HomVar variants respectively (with labels "", "Het", and "HomVar" respectively). 
Note also that SNP and INDEL statistics are automatically computed separately in the workflow, so there's no need to split 
them in your expressions. **WARNING**: Despite the examples in the `bcftools` documentation, you must follow the default 
formatting and use double quotes around the WDL string, with single quotes around internal expressions when appropriate.
* `bcf_labels` - (default = `[]`) a list of labels to match the input `bcf_selectors` list, to appear next to statistics 
computed on the corresponding subset in the output data tables.

* `evaluation_bed` - if provided, only variants inside this bed file will be counted in statistics. Note this is slightly 
different than using a `strat_interval` because `vcfeval` will process variants on the boundary slightly differently here 
(see the RTG manual for details). This is intended to be used when one has a "high confidence" region for known truth samples. 
**This must be a bed file!**
* `score_fields` - (default = `["GQ"]`) an optional list of fields to use for stratifying ROC data. For each entry, the 
`vcfeval` will be run to produce the proper ROC data with the given string input to `--vcf-score-field`. Example values 
include `"GQ"`, `INFO.<NAME>`, and `FORMAT.<NAME>`. **WARNING**: Some values might not make sense to stratify when 
interested in *both* precision and recall statistics. For example, with the value `FORMAT.DP` it might make sense to 
evaluate precision of calls made by depth in the callset, but often the baseline/truth doesn't have any intrinsic DP 
associated with each variant, and so any present in the baseline may lead to misleading recall statistics. See the 
section on [Interpreting Output Data](#interpreting-output-data) below for details.

* `experiment_label` - an optional string to populate "Experiment" column in final output tables. Useful for incorporating 
Terra table information about experimental groups to do grouped analyses downstream, or for easier integration into provided 
visualizer tools.
* `extra_column` - an optional string title for an additional column to add to final output tables, for downstream analysis.
* `extra_column_label` - values to populate in the optional `extra_column`, if specified. To make downstream analysis easier 
by incorporating Terra table information, e.g. `MEAN_COVERAGE`.

* `NULL_FILE` - do NOT overwrite. This input exists for technical reasons due to WDL lacking good handling for null File 
types in its current version.

In addition, there are some `vcfeval` options that can be toggled using the Boolean inputs `require_matching_genotypes`, 
`enable_ref_overlap`, and `passing_only`. Check the RTG manual for more details on how these can affect output. For most users, 
these can be safely left as default.

Each task also has inputs to control the amount of disk, CPU, and memory allocated. These can be toggled to optimize 
performance on different inputs.


### Interpreting Output Data

There are four files output by the workflow. They are:

* `simple_summary` - This file is approximately the same as the usual outputs from the old benchmarking script. It contains 
counts of the usual TP, FP, etc. statistics, along with computed Precision, Recall, and F1_Score stats.
* `combined_IDD` - (InDel Distribution) This file contains distributions of the TP, FP, etc. statistics by length of INDEL, 
where negative values correspond to deletions. For convenience, the column `INDEL_Type` labels the variant class as `Ins` or 
`Del` accordingly, so users can readily group and process their own statistics per class.
* `combined_ST` - (Substitution Types) This file contains distributions of the TP, FP, etc. statistics by SNP substitution 
(e.g. `A>C`, `A>G`, etc.). For convenience, the column `Substitution_Type` labels the variant class as `Ti` or `Tv` depending 
on if the substitution is a Transition or Transversion accordinly, so users can readily group and process their own 
statistics per class.
* `combined_ROC` - This file is a concatenation of all of the ROC outputs by `vcfeval`, with labels to separate out the 
stratifiers, variant types, etc. A `score_field` column records the user input score label for the given row.

We follow the same convention as RTG when computing Precision / Recall statistics, namely we use:

$$ Prec = \frac{TP_{call}}{(TP_{call} + FP)}, $$

$$ Recall = \frac{TP_{base}}{(TP_{base} + FN)}, $$

and

$$ F_1 = 2 \cdot \frac{Prec \cdot Recall}{Prec + Recall}.$$

#### Some Warnings on Interpreting Output Statistics

Although there is a lot of possible functionality afforded by the extensive input options above, one should be careful about 
interpreting the results in some cases. We provide a few examples here as warnings against turning a bunch of options on and 
focusing on statistics that may end up being misleading if not carefully analyzed, which is inherent to the benchmarking process.

##### 1. Subsetting regions can cause greater FPs/FNs.

Due to differences in representing equivalent variants, it's possible that subsetting can lead to false classification of 
variants as mistakes when they should be considered TPs. Suppose there is a variant $A$ in position $X$ in the call VCF, 
and a variant $B$ in position $X+Y$ in the base VCF. Normally, `vcfeval` will resolve the underlying haplotypes to compare 
if these produce equivalent sequences. However, if your input intervals mandate including only those variants at a position 
less/greater than somewhere between $X$ and $X+Y$, then only one of these variants will be visible during the `vcfeval` 
processing, and you won't get matching haplotypes (when normally you should).

To mitigate this, the `interval_padding` input is a dial that can be turned to see how variants near the boundaries of 
your interval files may be influencing the end statistics. Seeing large jumps when toggling this might be indicative of 
lots of action on the edges, and you might consider changing your interval files (or doing away with them altogether, and 
only using the `evaluation_bed` input which doesn't subset, but performs other advanced counting).

##### 2. Using ROC scores / variant selectors that are asymmetric can be misleading.

The inputs allowed for stratifying ROC scores and subsetting variant types with `bcftools` filtering selectors is extremely 
flexible, but this comes at the risk of producing statistics that may not make sense if not done carefully. Here's a good 
example to keep in mind when thinking about if your data might be susceptible to this.

Suppose you have a call VCF you'd like to benchmark against some "truth" base VCF. You want to understand how the depth 
(DP) of your callset correlates with better accuracy in the variants it calls. You put "DP" in as a `score_field`, and 
produce some ROC outputs. In this case, the Recall statistics output would not be useful, because "truth" VCFs are often 
produced using very large coverage. Because DP is not "intrinsic" to the variants themselves, but rather to the sequencing 
process, there's no canonical way to measure to what extent you accurately "recall" variants above a given DP threshold. 
The precision statistics, however, make sense to interpret above given DP thresholds in the call VCF because they are 
derived from TP and FP alone, which are labeled based on the call variants. In this way, you can track how accuracy 
(hopefully) improves when increasing DP in the calls, but should ignore the DP values in the base VCF.

On the other hand, consider another situation where it makes sense to subset by variant type, when the situation is 
sequencing-independent and more symmetric. Suppose you have a variant annotation that labels SNPs occuring inside 
reference homopolymers by the length of the homopolymer. Because this only depends on the reference, this is a much 
safer variant selector to subset on (maybe by setting some thresholds using the `bcftools` expression language) because 
variants occuring in either the call or base VCF should be on equal footing.

Note in the case that a ROC score field exists in the call VCF but not the base, you can still get ROC statistics output 
from `vcfeval`. See the RTG manual for details on how it interpolates values for the base along the different thresholds.

---

## Structure of the Workflow

The workflow has the following main high-level tasks:

1. **Subset** the pair of VCFs given to each interval file provided.
2. **Classify Variants** in both VCFs for each subsetted pair using `vcfeval`.
3. **Compute Statistics** on labeled variants in different subsets of types of variants.
4. **Combine Statistics** into a few final files for output.

Some details are provided for each step below.

### Subset

For each provided interval file, we use `gatk SelectVariants` to subset the provided VCFs. An optional padding amount 
can be provided to help control for a potential large number of variants on the boundaries. This always includes a special 
case where no subsetting is performed for "genone-wide" performance to be captured alongside the subsets provided. If a 
multisample VCF is provided, users must specify the sample name of the desired sample to proceed to compare in the next 
steps. Extraction of that sample is performed during this step.

### Classify Variants

For each subsetted pair from the last step, we run `vcfeval` to perform haplotype-aware variant comparison. This is a 
sophisticated way to comparing variants across VCFs which accounts for potentially different methods of representing 
haplotype-equivalent variants. Each variant in the files can be categorized as true positive (TP), false positive (FP), 
and false negative (FN) based on how the call variants compare with the base. The command outputs four VCF files containing 
each of these cateogies (treating TPs in call as possibly different from TPs in base) which then get passed to downstream 
analysis tools.

In addition, the command outputs ROC-style statistics split by user-input score fields (e.g. stratify TP, FP, FN counts 
by GQ scores, etc.). Users are allowed to specify a list of such score fields (e.g. GQ, QUAL, etc.) and the workflow 
will create ROC files for each. Some data cleaning is performed in this step on the ROC files output by `vcfeval`.

### Compute Statistics

In this step, we take the VCF subsets of TPs, FPs, and FNs and run `bcftools stats` on them to collect an array of 
statistics, particularly the count of SNPs and INDELs in each of these categories. This is performed once for each user-input 
variant selector, which must be communicated using the `bcftools` expressive filtering language. A large amount of data 
cleaning is also done in this step on the outputs to prepare for them to be compiled into the final outputs later.

### Combine Statistics

This step collects all of the statistics files generated previously over the different interval files and variant selector 
conditions and combines them into one large file. Each row will contain unique identifiers to explain the interval file 
subsetted to and any variant selector conditions used, along with user-provided input sample names, so that users can 
safely concatenate outputs across different runs of the workflow for multi-sample evaluation analysis. There is also a 
little bit more cleaning to create a "simple summary" output, which is reported alongside the more indepth statistics 
collected before.

---

## Comparison with BenchmarkVCFs

Although this tool was designed to be an improvement to the BenchmarkVCFs workflow, there are a few key differences that 
are important both for usage and functionality, which are briefly summarized below for users familiar with the old tool. 
In other words, this a quick guide for porting over inputs / workflows depending on BenchmarkVCFs to use this workflow instead.

### Output Differences
* The new `simple_summary` output is meant to closely mimic the `summary` output of BenchmarkVCFs. The key differences are:
    * The new file is a tsv while the old file is a csv.
    * The new file uses `TP_Call` in place of `TP_Eval` in the column names. In addition, the new file swaps `Call_Name` 
for `Name` and `Base_Name` for `Truth_Set`.
    * A few columns from the old file have been completely dropped. These are:
        * `Comparison_Engine` - the new tool only supports using `vcfeval` to compare VCFs.
        * `IndelLength` - users looking to break down statistics by INDEL length are encouraged to use the new `combined_IDD` output.
        * `IGV_Session` - to save some time, this feature of adding IGV session links has been dropped. This may be added back as an option in the future.
        * `Summary_Type` - this always output `summary` for a value, except in some INDEL stat cases. Because INDEL stratifying is handled by a new separate output, this column was dropped.
        * `UNK` - this previously counted variants outside of your evaluation region. Because this can be inferred from the other statistics reported, this was dropped to save time.
    * Some new default variant categories are included in the new summary, whose counts are computed for free by `bcftools stats`. 
Note when a site has multiple variants (e.g. a SNP and an INDEL), it increments both counters. These are:
        * MNPs
        * Other
        * MA (multiallelic sites)
        * MASNP (multiallelic SNP sites)
    * Subsetting by variant type is now handled by `bcf_selectors` using the `bcftools` filtering expression language, as opposed to GATK jexl expressions (which also did not work properly in the previous WDL in many circumstances).
    * Subsetting by interval files no longer automatically includes subsetting on their complements. Because many interval 
  files are specialized, there is major overlap in their complements and often not as much interpretability. Instead, users are 
  encouraged to create complements of their interval files using any tools available when this data is desired, and then include 
  them on equal footing with the other interval input files.
    * The default `bcf_selectors` automatically result in breakdowns by all vs Het vs Hom sites. The data processing is 
  handled so the `bcf_labels` are concatenated in front of the SNP / INDEL / etc labels, so you still get Type categories like 
  "HetSNP", etc. In particular, this means e.g. that the old `HetIndel` becomes `HetINDEL` in the new version.
    * If you run the old and new tools on the same inputs, you might notice some very small changes in the output statistics. 
  There are few factors here. One is that the old workflow is using an older version of RTG, and some changes to `vcfeval` 
  have been made since then which get reflected in slightly different stats. Another, more important difference is that 
  `bcftools` is used to generate a number of the statistics, which may calculate variants slightly differently than the 
  old tools. (E.g. Try counting the number of SNPs in a VCF using `gatk`, `bcftools`, and `rtg`, and you'll often get slightly 
  different numbers.) If you notice a large divergence in the final statistics, in particular enough to cause differences in scientific conclusions, please open an [issue](https://github.com/broadinstitute/palantir-workflows/issues).
* There are three new file outputs which help give granular statistics, taking into account INDEL length distribution, types of SNPs, etc. You can learn more about these in the [Usage](#usage) section.

### Internal Differences
* The number of tasks run in a typical run of the new workflow should be much lower than the old. This should translate to a significant decrease in 
runtimes, though that can vary with cloud weather.
* The old `CrosscheckFingerprints` task (`MatchEvalTruth`) in the BenchmarkVCFs workflow has been removed for now. Any good 
benchmarking workflow should include a fingerprint check at some step, but it wasn't clear if that should be included in *this* WDL, 
or a wrapper analogous to FindSamplesAndBenchmark, which has to run fingerprinting to pair up samples regardless. It is likely 
this will be added back as an (optional, but highly recommended) step in this WDL in the future.