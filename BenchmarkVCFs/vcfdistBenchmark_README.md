# vcfdistBenchmark Guide

This guide describes an overview of the [vcfdistBenchmark](vcfdistBenchmark.wdl) WDL used for evaluating VCF performance 
against a known baseline.

Some useful resources referenced below:

---

## Acknowledgements

This tool was created to run vcfdist, created by Tim Dunn, on Terra.

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
To run this file, upload it to a workspace on Terra, and run it on a table with the input information it needs. More information below,
but please note that **both input vcfs need to be phased**. As of 2023-08-24, there are **no warnings about data inputs that could be unphased**.
vcfdist will not crash, but your accuracy will likely be affected. Also, vcfdist ignores (as of 2023-08-24) whether the variants' haplotypes contain
a ``/`` or a ``|`` to improve compatibility with tools that have non-standard outputs.

### Inputs

The following are the top-level inputs for the workflow. We do not use the `vcfeval` terminology, which refers to the paradigm 
of a "base"line VCF and a "call"set VCF benchmarked against it. Instead, Tim Dunn's terminology is used. The "base"line VCF is here called
the "truth" VCF, and the "call"set VCF benchmarked against it is referred to as the "query" or "eval" VCF.

#### Necessary Inputs
* `truth_vcf` - the VCF to use as baseline. **Must be locally phased!**
* `eval_vcf` - the VCF ("callset") to be evaluated. **Must be locally phased!**
* `bed_file` - only variants inside this bed file will be counted in statistics. **This must be a bed file.**
* `fasta_file` - reference fasta file matching the input VCFs.

#### Optional Inputs
* `cpu` - optional, integer input. Default is `2`
* `disk_size_gb` - optional, integer input. Default is the size of the truth VCF in GiB, plus 10, rounded up.
* `mem_gb` - optional, integer input. Default is `10`.
* `preemptible` - optional, integer input. Default is `1`.
* `docker` - optional, string input. Default is `us.gcr.io/broad-dsde-methods/vcfdist:v0.1`.

### Interpreting Output Data

There are five files output by the workflow. They are:

* `prs_tsv` - *"High-level precision/recall overview of SNP/INDEL performance. For each category (SNP/INDEL), there is one line for performance at the chosen minimum quality score, and one line for the quality score threshold that results in the best performance."* -vcfdist documentation
* `precrec_tsv` - *"For each category (SNP/INDEL), there is one line reporting the precison/recall performance at each possible quality score."* -vcfdist documentation
* `query_tsv` - *"Reports detailed information regarding each variant"* (-vcfdist documentation) for the query vcf file.
* `truth_tsv` - *"Reports detailed information regarding each variant"* (-vcfdist documentation) for the truth vcf file.
* `vcfdistsummary` - A vcf summarizing the variants, output by vcfdist.

#### Some Warnings on Interpreting Output Statistics

##### 1. Because data needs to be phased, the quality of the benchmarking will be affected by the sequencing quality, the 
read quality, and the phasing.

##### 2. [].
[]

##### 3. [].

[]

---

## Structure of the Workflow

The workflow has one high-level task:

1. **Subset** the pair of VCFs given to each interval file provided.
2. **Classify Variants** in both VCFs for each subsetted pair using `vcfeval`.
3. **Compute Statistics** on labeled variants in different subsets of types of variants.
4. **Combine Statistics** into a few final files for output.


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