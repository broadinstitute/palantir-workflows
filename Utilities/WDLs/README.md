# Utility WDLs

This directory contains a collection of miscellaneous WDLs useful for some small tasks. Check below for documentation
on each.

* [CollectBenchmarkSucceeded](#collectbenchmarksucceeded)
* [Dipcall](#dipcall)
* [IndexCramOrBam](#indexcramorbam)

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

