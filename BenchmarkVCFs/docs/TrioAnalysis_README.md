# TrioAnalysis Guide

## Overview
The `TrioAnalysis` workflow is designed to analyze the genetic data of a child and their parents to compute Mendelian concordance statistics. The workflow takes VCF files for the child, father, and mother, along with reference genome files, and optionally annotates the VCF with regions from BED files. The final output includes annotated VCF files and summary tables of Mendelian violations and uncertainties.

## Inputs
- **`File child_vcf`**: VCF file for the child.
- **`File child_vcf_index`**: Index file for the child's VCF.
- **`String sex_of_child`**: Sex of the child, either "MALE" or "FEMALE".
- **`File father_vcf`**: VCF file for the father.
- **`File father_vcf_index`**: Index file for the father's VCF.
- **`File mother_vcf`**: VCF file for the mother.
- **`File mother_vcf_index`**: Index file for the mother's VCF.
- **`File ref_fasta`**: Reference genome FASTA file.
- **`File ref_index`**: Index file for the reference genome.
- **`Array[File]? bed_files`**: Optional array of BED files for annotation.
- **`Array[String]? bed_labels`**: Optional array of labels corresponding to the BED files.

The `MakeSummaryTable` task also takes a `format_value` input for a FORMAT field from the VCF to create a CDF of mendelian violation counts from. This is optional and defaults to "GQ". The `format_agg` value sets how to combine these values across the trio and can take the values: `min` (default), `max`, `mean`, `median`, `sum`. So for example, the default behavior would be to make a cumulative distribution of counts of mendelian violations stratified by the minimum GQ value in the trio for each site. One would expect with accurate sequencing/data processing, this count should approach the de novo rate as the threshold on min GQ increases. 

## Outputs
There are two types of variant triples of concern the `rtg mendelian` tool will flag: "mendelian violations" and "mendelian uncertainties."

* **Mendelian violation: a variant triple inconsistent with standard inheritance. E.g. father: 0/0, mother: 0/1, child: 1/1.**
* **Mendelian uncertainty: a variant triple that could be a mendelian violation, but missing some genotypes to know for sure. E.g. father: ./., mother: 0/1, child: 1/1; or father: ./., mother: 1/1, child: 1/1.**

Because of this, you may get more useful counts of mendelian uncertainties using a gVCF (perhaps with a GQ threshold) than a standard VCF. These two types of variant triples are split into the following outputs.

- **`File mendelian_output_vcf`**: The final VCF file after Mendelian concordance analysis and optional annotation.
- **`File mendelian_output_vcf_index`**: Index file for the final VCF.
- **`File? mendelian_violation_table`**: Table of Mendelian violations.
- **`File? mendelian_uncertainty_table`**: Table of Mendelian uncertainties.
- **`File? mendelian_violation_counts`**: Counts of Mendelian violations across provided regions.
- **`File? mendelian_uncertainty_counts`**: Counts of Mendelian uncertainties across provided regions.
- **`File? mendelian_violation_cdf`**: Cumulative distribution function (CDF) of Mendelian violations for each provided region.
- **`File? mendelian_uncertainty_cdf`**: Cumulative distribution function (CDF) of Mendelian uncertainties for each provided region.

## Workflow Description
1. **ComputeMendelianConcordance**: This task normalizes and merges the VCF files for the child, father, and mother, and computes Mendelian concordance statistics using the `rtg mendelian` tool. The output is a VCF file with Mendelian concordance annotations.
2. **AnnotateVcfRegions** (optional): If BED files and labels are provided, this task annotates the VCF file with regions from the BED files using `bcftools annotate`.
3. **MakeSummaryTable**: This task generates summary tables of Mendelian violations and uncertainties, including counts and cumulative distribution functions (CDFs).
