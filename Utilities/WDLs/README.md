# Utility WDLs

This directory contains a collection of miscellaneous WDLs useful for some small tasks. Check below for documentation on each.

* [CollectBenchmarkSucceeded](#collectbenchmarksucceeded)
* [Dipcall](#dipcall)
* [IndexCramOrBam](#indexcramorbam)
* [DownsampleAndCollectCoverage](#downsampleandcollectcoverage)

## CollectBenchmarkSucceeded

### Summary

When running the `FindSamplesAndBenchmark` wdl with many samples, it sometimes happens that a few fail in Terra while most succeed.
Unfortunately, this means that the outputs of benchmarking the successful ones don't get compiled into one convenient .csv file to use for data analysis.
If you don't mind sacrificing the few that failed, or want to get started analyzing the successful ones ASAP, this wdl will automatically collect
the successful outputs and aggregate them into one .csv, similar to the last task of the benchmarking wdl.

### Inputs

* `namespace`: the first personalized part of your workspace URL; e.g. if you see `<my_project>/<my_workspace>` at the top
  in Terra, then this should be `<my_project>` as a string.
* `workspace_name`: specific name for your workspace, e.g. `<my_workspace>` in the last example.
* `submission_id`: the submission id for the `FindSamplesAndBenchmark` run, found from the "Job History" tab.

## Dipcall

### Summary

This WDL is a modified version of the Dockstore version of the [Dipcall](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/QC/wdl/tasks/dipcall.wdl)
pipeline. This workflow takes in a diploid assembly and calls variants, creating a VCF against your chosen reference.
This modified version allows for you to specify custom PAR regions for your reference so you can call haploid variants
when appropriate. Some data cleaning and indexing of the output VCF is also performed.

### Inputs

* `assemblyFastaPat`: the haploid assembly fasta for paternally inherited chromosomes
* `assemblyFastaMat`: the haploid assembly fasta for maternally inherited chromosomes
* `referenceFasta`: the reference to call variants against
* `isMaleSample`: set true if you would like to make haploid calls on X/Y outside the PAR region
* `custom_PAR_bed`: bed file denoting pseudoautosomal (PAR) regions for your reference
* `sample_name`: name to put for your sample in final output VCF
* `referenceIsHS38`: set true (default) if using hg38 reference


## IndexCramOrBam

### Summary 

Use this WDL to index a CRAM or BAM file, using `samtools`. The type is inferred using the file extension (either `.cram` or `.bam`). 


## DownsampleAndCollectCoverage

### Summary

The idea of this WDL is to do everything you need for a standard downsampling experiment. It takes in either CRAM or BAM files and downsamples them either according to a defined downsampling ratio or to a desired target coverage. If no downsampling ratio is defined then it will run `ColectWgsMetrics` to get the original mean coverage and determine the downsampling ratio based on that coverage and the desired target coverage. After downsampling the workflow will run CollectWgsMetrics once more and output the mean coverage of the downsampled CRAM file. This provides feedback with respect to the target coverage, because downsampling is always associated with some uncertainty.

### Inputs 
* `File input_cram`: Input BAM or CRAM
* `File input_cram_index`: Index file for input BAM or CRAM
* `File ref_fasta`: Reference FASTA
* `File ref_fasta_index`: Reference FASTA index
* `Float? downsample_probability`: Downsampling ratio. If not provided, the ratio will be determined based on the `target_coverage`.
* `Float? target_coverage`: Target mean coverage for the downsampled CRAM file. **In order to use this input, do not provide `downsample_probability`, otherwise, that value will be used for downsampling.** 
* `File? coverage_intervals`: If provided, the output downsampled mean coverage will be calculated based on these intervals. Additionally, these intervals will be used to calculate the original coverage if `target_coverage` is used.
* `String downsample_strategy = "ConstantMemory"`: See [DownsampleSam documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832708637467-DownsampleSam-Picard-).
* `Int read_length = 150`: See [CollectWgsMetrics documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832707851035-CollectWgsMetrics-Picard-).
* `Boolean use_fast_algorithm = true`: See [CollectWgsMetrics documentation](https://gatk.broadinstitute.org/hc/en-us/articles/13832707851035-CollectWgsMetrics-Picard-).
* `String docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"`: Docker to use for both CollectWgsMetrics and DownsampleSam
* `File? picard_jar_override`: If provided, a Picard JAR file to use for both CollectWgsMetrics and DownsampleSam instead of the `gatk` command.
* `Int preemptible = 1

### Outputs
`File downsampled_cram`: Downsampled CRAM
`File downsampled_cram_index`: Downsampled CRAM index
`Float downsampled_mean_coverage`: Mean coverage over the `coverage_intervals` (or the whole genome if not provided) for the downsampled CRAM file
`Float? original_mean_coverage`: The original mean coverage over the `coverage_intervals` (or the whole genome if not provided) of the input CRAM file, if `target_coverage` was used
