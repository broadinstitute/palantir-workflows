# Utility WDLs

This directory contains a collection of miscellaneous WDLs useful for some small tasks. Check below for documentation on each.

* [CollectBenchmarkSucceeded](#collectbenchmarksucceeded)
* [CreateIGVSession](#createigvsession)
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

## CreateIGVSession

### Summary

This workflow takes in optional lists of BAMs, VCFs, and interval files (`.interval_list` or `.bed`) and combines them together
into an IGV session .xml file. A reference must be provided, either by a hardcoded string ("hg38" or "hg19"), or by providing a 
path to the desired fasta. Input files are interpreted as WDL strings, so no localization occurs. Bucket paths are output in the .xml
session, so IGV will stream them directly from the cloud. This task is useful to add to the end of workflows that output lots of files
you might want to visualize together for analysis or debugging.

### Inputs

* `bams`: (optional) list of BAMs/CRAMs to add to session.
* `vcfs`: (optional) list of VCFs to add to session.
* `interval_files`: (optional) list of `.interval_list` or `.bed` files to add to session.
* `reference`: reference to use in IGV; must be either a `.fasta` file or one of the values: "hg38" or "hg19".
* `output_name`: (default = "igv_session") name for the output .xml file.

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


## IntervalList2Bed

### Summary

This WDL takes in a list of interval files (either `.bed` or `.interval_list`) and converts the `.interval_list` files into `.bed`. The WDL checks if any of the provided files has a `.interval_list` extension, and then will call a conversion task on it if so. This means if all the files provided are `.bed`, then no tasks will be called, and the original list will be returned. This allows you to drop this task in to your workflows to extend pipeline functionality from accepting `.bed` inputs to also handle `.interval_list` files without penalizing users who provided `.bed` files with unnecessary extra tasks, which is ideal as many tools require specifically `.bed` lists.

If labels are provided, they will be returned in the new order of the `bed_files` output, which may be different than the originally given order. If labels are not provided, a list of `basename`s for the input files will be returned, in the correct order.

### Inputs
* `interval_files`: a list of `.bed` and/or `.interval_list` files
* `interval_labels`: an optional list of string labels to use for the corresponding interval file

### Outputs
* `bed_files`: a list of `.bed` files converted to the given inputs; note the order may have changed from the given list
* `bed_labels`: a list of labels for the `.bed` files corresponding to the user inputs, or the basename of the input files if the user did not provide any labels. Note the order may have changed, but the position of a label corresponds to the position of the file in `bed_files`.


## DownsampleAndCollectCoverage

### Summary

The idea of this WDL is to do everything you need for a standard downsampling experiment. It takes in either CRAM or BAM files and downsamples them either according to a defined downsampling ratio or to a desired target coverage. If no downsampling ratio is defined then it will run `ColectWgsMetrics` to get the original mean coverage and determine the downsampling ratio based on that coverage and the desired target coverage. After downsampling using `DownsampleSam` the workflow will run `CollectWgsMetrics` once more and output the mean coverage of the downsampled CRAM (or BAM) file. This provides feedback with respect to the target coverage, because downsampling is always associated with some uncertainty. If `fail_if_below_coverage` is set, the workflow will fail if that downsampled mean coverage is below the provided threshold.

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
