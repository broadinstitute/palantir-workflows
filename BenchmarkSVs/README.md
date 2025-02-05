# BenchmarkSVs

This document describes the usage and outputs of the BenchmarkSVs WDL along with an accompanying plotting script for creating a simple dashboard for visualizing the outputs, called the SVisualizer.

## Table of Contents

1. [BenchmarkSVs WDL Overview](#benchmarksvs-wdl-overview)
   1. [Summary](#summary)
   2. [Inputs](#inputs)
   3. [Outputs](#outputs)
2. [SVisualizer](#svisualizer)
   1. [Dependencies](#dependencies)
   2. [Usage](#usage)
   3. [Optional Data Gathering Script](#optional-data-gathering-script)
3. [The CleanSVs WDL](#the-cleansvs-wdl)
4. [A Note on Wittyer](#a-note-on-wittyer)

## BenchmarkSVs WDL Overview

### Summary

There are two primary workflows contained in evaluating performance of SV VCFs. Both call some common tasks from the [tasks](tasks.wdl) file. The main workflows are:
1. [QcSVs.wdl](QcSVs.wdl): Collect QC Metrics -- uses `bcftools` to generate some statistics on a cohort VCF, like AF, HWE p-values, etc. per site.
2. [BenchmarkSVs.wdl](BenchmarkSVs.wdl): Run Truvari -- uses `truvari` to compare your input `comp_vcf` against the `base_vcf` provided. 

In the latter, statistics are collected over various subsets of your data, including relative to a list of interval files the user can supply.

The final outputs include tables of cleaned data that can then be used for your downstream analysis, either on your own or using the provided [SVisualizer](#svisualizer).

### BenchmarkSVs Workflow

#### Inputs

- `base_vcf` - proxy for "truth" VCF to use; must contain `END`, `SVTYPE`, and `SVLEN` INFO fields.
- `base_vcf_index` - index for `base_vcf`
- `base_sample_name` - sample name to use when comparing via a cohort VCF
- `comp_vcf` - the main input VCF you want to QC/benchmark; must contain `END`, `SVTYPE`, and `SVLEN` INFO fields.
- `comp_vcf_index` - the index for `comp_vcf`
- `comp_sample_name` - sample name to use when comparing via cohort VCF
- `experiment` - (optional) label to populate output tables with. Useful for running multiple iterations of the workflow, and then concatenating outputs to use this label for downstream clustering in plots.
- `ref_fasta` - reference fasta file
- `ref_fai` - reference fasta index file
- `evaluation_bed` - (optional) a bed to subset both inputs VCFs to before performing any of the following analysis
- `evaluation_pct` - (default 1) a float between 0 and 1 for the percent overlap with `evaluation_bed` required when subsetting to keep
- `bed_regions` - (optional) bed files to collect stats relative to
- `bed_labels` - (optional) labels to use in output tables for corresponding `bed_regions` files
- `breakpoint_padding` - (default = 20) amount of bases to add to the breakpoint edges of SVs for computing overlap stats
- `svlen_bin_cutoffs` - (default = [100, 250, 1000, 2500, 10000, 25000]) bin edges to use for categorizing SVs by length
- `create_igv_session` - (default = true) toggle whether to run [CreateIGVSession](/Utilities/WDLs/CreateIGVSession.wdl)
- `optional_igv_bams` - (optional) input for [CreateIGVSession](/Utilities/WDLs/CreateIGVSession.wdl)

In addition to these, specific tasks expose flags for certain tools. The `RunTruvari` task, for example, has inputs for toggling all the available flags for the tool. See the [documentation](https://github.com/acenglish/truvari/wiki/bench) for that tool for more details. 

#### Outputs

- `truvari_bench_summary` - table containing summary stats produced by Truvari, combined over interval files provided
- `truvari_fn_closest` - table containing data on the FNs marked by Truvari along with their (default = 3) closest variants in the opposing file
- `truvari_fp_closest` - table containing data on the FPs marked by Truvari along with their (default = 3) closest variants in the opposing file
- `igv_session` - output from [CreateIGVSession](/Utilities/WDLs/CreateIGVSession.wdl) if toggled on

### QcSVs Workflow

#### Inputs

- `vcf`: input file for SV calls
- `vcf_index`: index for `vcf`

For documentation on other arguments, see the above on `BenchmarkSVs`.

#### Outputs

- `qc_summary` - table containing cohort level stats per site, including AF, HWE, etc.

## SVisualizer

The SVisualizer is a python script based on [Dash](https://plotly.com/dash/) and [Quickboard](https://github.com/broadinstitute/quickboard/) to create an interactive dashboard using the outputs of the above WDL. Currently, the board assumes file outputs will exist for each of the three main components of the workflow. To use, follow the instructions below.

### Dependencies

Follow these steps to ensure your environment is prepared.

1. Make sure you have the Python package `quickboard` installed. It's recommended to pin the version to avoid any breaking changes in the future. The version `quickboard==0.4.0` should work with this current version.
2. Copy all the files from the [SVisualizer](SVisualizer) directory to a fresh directory on your local system. 
3. Copy all the files from the WDL outputs into a subdirectory called `wdl_outputs`. You can also use the gathering script [provided](#optional-data-gathering-script) to make this step easier.

### Usage

Once you've followed all the steps from the [Dependencies](#dependencies), using the script is straightforward. You can open and run the `.ipynb` file provided in Jupyter. This will start a local Flask server (hosted on port 8050 by default) you can navigate to in order to interact with the dashboard. It will be automatically populated using the WDL output data in the `wdl_outputs` directory. The notebook version has some markdown headers to make it easier to navigate the different blocks of code, if you wish to modify it for your own purposes. It also makes it easier to see what variables/settings can be toggled at the top of the page.

### Optional Data Gathering Script (Legacy)

In the legacy version of a single workflow, there are a lot of files output, and you might want to use the provided `gather_terra_data.py` script to automatically set up your environment. This works using the Firecloud API to grab your files from a Terra submission and put them in the right location, combining outputs from separate workflows run under the same submission. To use it, you need to:

1. Make sure you have the Firecloud API Python package installed, e.g. `pip install firecloud`.
2. Make sure you have `gsutil` installed and configured with the correct permissions relative to your Terra workspace.
3. Run `python gather_outputs.py <TERRA_NAMESPACE> <TERRA_WORKSPACE_NAME> <SUBMISSION_ID>`.

This will create the `wdl_outputs` directory with the right files inside. You can find the command arguments by navigating to the Job History tab and clicking on your desired run. You can use e.g. the Terra URL to figure out most of the arguments. For example, the URL should look like:
```
app.terra.bio/#workspaces/<TERRA_NAMESPACE>/<TERRA_WORKSPACE_NAME>/job_history/<SUBMISSION_ID>
```

The script works by downloading the combined files `.tar` directory output from the workflow, extracting it to a local temp directory, and then combining these outputs across all successful workflow runs under the same submission. Note separate "Experiment" label inputs can be used to distinguish parts of the output tables corresponding to different runs. The script naively concatenates the files across the runs, so doing this for a large number of workflows inside the same submission could result in large tables, which may not be ideal for the dynamic SVisualizer.


## The CleanSVs WDL

Some tools are particular about the format required to process SV VCFs, like Wittyer (deprecated). As a helper script, the CleanSVs WDL is provided to reformat SV VCFs to make them more amenable for certain types of analysis. Its functionality and usage are briefly described here.

### Inputs

- `input_vcf` - the VCF to clean
- `input_vcf_index` - (optional) the index for the input VCF
- `output_name` - name for the output VCF file
- `min_size` - (default = 50) minimum (absolute value) size of variants to keep
- `split_MA` - (default = true) split multiallelic sites into biallelic sites
- `add_annotations` - (default = true) add annotations like `END`, etc. required by the benchmarking pipeline if not already included
- `normalize_and_clean` - (default = true) remove `SVTYPE`s outside specified list, and convert missing filter to `PASS` optionally; also includes option to set min length bound and replace all filters with `PASS` in task inputs
- `convert_to_abstract` - (default = true) convert SVs represented by long sequences into variants with symbolic alleles (`INS` or `DEL` only)

### Outputs

- `output_vcf` - final cleaned VCF
- `output_vcf_index` - index for final output

## A Note on Wittyer

Previous iterations of this WDL used the Wittyer tool for benchmarking, however this functionality has been deprecated in favor of using Truvari to perform the comparisons for greater flexibility and less redundancy. However, legacy code is provided both in the [WittyerTasks.wdl](WittyerTasks.wdl) and [wittyer_tabs.py](SVisualizer/wittyer_tabs.py) files. These are less polished and tested than the rest of the code, but may prove to be useful for anyone looking to experiment with this benchmarking engine in the future.