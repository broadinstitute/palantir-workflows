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

## Comparison with BenchmarkVCFs (and the tool vcfeval): Some Warnings on Interpreting Output Statistics
* Since data needs to be phased, it's important to remember that the quality of the benchmarking will be affected 
by the sequencing quality, the read quality, and the phasing.
* The calculated precision and recall values output in `prs_tsv` and `precrec_tsv` includes the effects of partial credit. 
To get the precision and recall statistics without partial matching, calculate them using the data in `precrec_tsv`
that corresponds to the quality score threshold that results in the best performance (the row with the larger MIN_QUAL); add the number of partial positives
shown for the truth vcf (TRUTH_PP) to the reported value of False Negatives (TRUTH_FN) to get the updated value of False 
Negatives (TRUTH_FN), and add the number of partial positives shown for the query vcf (QUERY_PP) to the reported value of 
False Positives (QUERY_FP) to to get the updated value of False Positives (QUERY_FP).
* **Because variants are outputted in biallelic rows, that affects the event counts (as each row is one event), and therefore the precision and recall values**
* **Because some INDELs in repetitive regions are replaced with SNPs and altered surrounding INDELs, that affects the event counts, and therefore the precision and recall values**
* **Please ensure that your data is locally phased!** Global phasing (as of 2023-08-24) is not required.
---

## Structure of the Workflow

The workflow has one high-level task: it runs vcfdist on the inputs and returns the output files listed above.



### Output File Structure Differences


---

### Internal Differences
* The number of tasks run in a typical run of the new workflow should be much lower than the old. This should translate to a significant decrease in 
runtimes, though that can vary with cloud weather.
* The old `CrosscheckFingerprints` task (`MatchEvalTruth`) in the BenchmarkVCFs workflow has been removed for now. Any good 
benchmarking workflow should include a fingerprint check at some step, but it wasn't clear if that should be included in *this* WDL, 
or a wrapper analogous to FindSamplesAndBenchmark, which has to run fingerprinting to pair up samples regardless. It is likely 
this will be added back as an (optional, but highly recommended) step in this WDL in the future.

---