# BenchmarkVCFs
**BenchmarkVCFs.wdl** produces a table of benchmarking statistics.  A skeleton (non-working) input json is included as **test.json**.  The following sets of arrays are linked (i.e. the first element in `Benchmark.stratIntervals` must correspond to the first element in `Benchmark.stratLabels`):
- Stratifier Intervals (optional): `Benchmark.stratIntervals`, `Benchmark.stratLabels`
- Variant Types (optional): `Benchmark.jexlVariantSelectors`, `Benchmark.variantSelectorLabels` 

There is no daisychaining of stratifiers, if you want to combine multiple stratifiers, you must do so youself and provide the combined interval list as a stratifier.  Confidence intervals and stratifiers may be specified either as Picard interval lists (ending in ".interval_list") or bed files (ending in ".bed").  Both file types can be included in the same run, and even in the same input array.  If the optional stratifier arrays are not included, the corresponding columns in the output tables will be populated with NA.  Variant type selection is performed using [GATK SelectVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php).  The JEXL indicated by `Benchmark.jexlVariantSelectors` will be passed to the `-select` option of SelectVariants.  For example, `vc.isIndel() && vc.getHetCount() == 1` will calculate metrics for all het indels.  Variant types are crossed with stratifiers, when both are provided.  Metrics are always caluclated for the basic SNP and INDEL variant types, along with Heterozygous and Homozygous variant SNPs and INDELs, regardless of whether specific variant type selections are given.  *Benchmark.referenceVersion* must be either "hg19" or "hg38", and correspond to the reference version of the fasta given in Benchmark.reference (this is needed for the igv xml).  By default stratification will also be performed over INDEL length.  This can be turned off with `Benchmark.doIndelLengthStratification`

The table of results will be found in the **call-CombineSummaries** bucket of the workflow as *summary.csv*.

# Output Format
The output csv file will have the following columns:
- **Name** : The name of the of the vcf being evaluated, as set in *Benchmark.evalLabels*
- **Truth_Set** : The name of the truth set being used, as set in *Benchmark.truthLabels*
- **Comparison_Engine** : The comparison engine used.  Current possibilities are GATK_GC for GATK/Picard GenotypeConcordance, Happy for hap.py, and VcfEval for vcfeval.
- **Stratifier** : The name of the stratifier used, if any, as set in Benchmark.stratLabels.  If no stratifier was used for this row value, will be *NA*
- **Summary_Type** : The type of summary this row refers to (with respect to indel length stratification).  Currently there are five options:
	1. *summary*: no stratification by indel length 
	2. *indel_fine*: indels have been binned by length with a bin width of 1 base
	3. *indel_coarse*: indels have been binned by length with a bin width of 5 bases
	4. *insertion* : all insertions
	5. *deletion* : all deletions
- **IndelLength** : The center of the indel length bin. If **Summary_Type** is *summary*, *insertion*, or *deletion* this value will be *NA*.  In order to know what indel lengths are included in this row, you must know both this value and **Summary_Type**.  For example, **IndelLength** 10 with **Summary_Type** *indel_fine* represents all indels with length 10, while **IndelLength** 10 with **Summary_Type** *indel_coarse* represents all indels with length 8,9,10,11, or 12.  Note that insertions are given positive **IndelLength** while deletions are given negative **IndelLength**.
- **Type** : If **Summary_Type** is *summary*, then this column tells whether the row if for SNPs, INDELs, or another more complicated type specified in `Benchmark.jexlVariantSelectors`.  If **Summary_Type** is not *summary*, the value will be *NA*.
- **TP_Eval** : Number of true positives in the vcf being evaluated.
- **TP_Base** : Number of true positives in the truth set.  Note that **TP_Eval** and **TP_Base** can be different due to differences in variant representation.  However, if the difference between them is large this is likely a sign of some problem.
- **FP** : Number of false positives (counted in the vcf being evaluated).
- **FN** : Number of false negatives (counted in the truth set).
- **UNK** : Number of variants outside the truth set high confidence region. At the moment this is defined only when **Comparison_Engine** is *Happy* and is *NA* for all others.
- **Recall** : **TP_Base**/(**TP_Base** + **FN**)
- **Precision** : **TP_Eval**/(**TP_Eval** + **FP**)
- **F1_Score** : Harmonic average of precision and Recall = 2(**Precision** x **Recall** )/(**Precision** + **Recall**)
- **IGV_Session** : The Google Bucket location of XML IGV session corresponding to this row.      
