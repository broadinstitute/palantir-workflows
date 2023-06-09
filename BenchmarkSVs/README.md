# BenchmarkSVs

This document describes the usage and outputs of the BenchmarkSVs WDL along with an accompanying plotting script for 
creating a simple dashboard for visualizing the outputs, called the SVisualizer. This is currently an experimental 
pipeline subject to lots of (potentially breaking) changes in the near future. 

## Table of Contents

1. [BenchmarkSVs WDL Overview](#benchmarksvs-wdl-overview)
   1. [Summary](#summary)
   2. [Inputs](#inputs)
   3. [Outputs](#outputs)
2. [SVisualizer](#svisualizer)
   1. [Dependencies](#dependencies)
   2. [Usage](#usage)
   3. [Optional Data Gathering Script](#optional-data-gathering-script)

## BenchmarkSVs WDL Overview

### Summary

This workflow allows you to collect some QC metrics on an input VCF and/or benchmark compare two SV VCFs. There are 
currently three key overall steps to the process, each optional depending on user inputs:
1. Collect QC Metrics -- uses `bcftools` to generate some statistics on a cohort VCF, like AF, HWE p-values, etc. per site.
2. Run Truvari -- uses `truvari` to compare your input `comp_vcf` against the `base_vcf` provided. 
3. Run Wittyer -- uses `witty.er` to compare your input `comp_vcf` against the `base_vcf` provided.

In the latter two steps, statistics are collected over various subsets of your data, including relative to a list of interval
files the user can supply.

The final outputs include tables of cleaned data that can then be used for your downstream analysis, either on your own
or using the provided [SVisualizer](#svisualizer).

### Inputs

- `base_vcf` - (optional) proxy for "truth" VCF to use
- `base_vcf_index` - (optional) index for `base_vcf`
- `base_sample_names` - (optional) list of sample names to use when comparing via a cohort VCF. Note the order must match
the order of `comp_sample_names` to get the correct matches. If not provided, only the first sample will be used.
- `comp_vcf` - the main input VCF you want to QC/benchmark
- `comp_vcf_index` - the index for `comp_vcf`
- `comp_sample_names` - (optional) list of sample names to use when comparing via cohort VCF. Note the order must match
  the order of `base_sample_names` to get the correct matches. If not provided, only the first sample will be used.
- `experiment` - (optional) label to populate output tables with. Useful for running multiple iterations of the workflow,
and then concatenating outputs to use this label for downstream clustering in plots.
- `perform_qc` - (default = true) decides whether to run QC task
- `run_wittyer` - (default = true) decides whether to run witty.er; note that witty.er expects certain types of SVs only
to be present, so make sure your VCF is compliant before running through this task. The types supported by witty.er are:
DEL, INS, DUP, INV, IntraChrBND, TranslBND. 
- `bed_regions` - (optional) bed files to collect stats relative to
- `bed_labels` - (optional) labels to use in output tables for corresponding `bed_regions` files
- `sv_types` - (optional) types of SVs to use for running Truvari comparisons across. Note an "ALL" type is automatically
computed in the workflow.

In addition to these, specific tasks expose flags for certain tools. The `RunTruvari` task, for example, has inputs for 
toggling all the available flags for the tool. See the [documentation](https://github.com/acenglish/truvari/wiki/bench) 
for that tool for more details. 

### Outputs

- `qc_summary` - table containing cohort level stats per site, including AF, HWE, etc.
- `truvari_bench_summary` - table containing summary stats produced by Truvari, combined over interval files provided.
- `truvari_fn_closest` - table containing data on the FNs marked by Truvari along with their (default = 3) closest variants 
in the opposing file.
- `truvari_fp_closest` - table containing data on the FPs marked by Truvari along with their (default = 3) closest variants
  in the opposing file.
- `truvari_fn_intervals` - table containing data on the FNs marked by Truvari along with information on their overlapping stats
with provided bed files
- `truvari_fp_intervals` - table containing data on the FPs marked by Truvari along with information on their overlapping stats
  with provided bed files
- `truvari_tp_base_intervals` - table containing data on the TPs in base marked by Truvari along with information on their 
overlapping stats with provided bed files
- `truvari_tp_base_intervals` - table containing data on the TPs in comp marked by Truvari along with information on their 
overlapping stats with provided bed files
- `wittyer_basic_stats` - table of basic stats produced by witty.er, including breakdowns over SV type, lengths, etc.
- `wittyer_truth_intervals` - table of variants in base with various labels annotated by witty.er
- `wittyer_query_intervals` - table of variants in comp with various labels annotated by witty.er
- `wittyer_nogt_intervals` - table of variants that were reported without GT for both base and comp according to witty.er


## SVisualizer

The SVisualizer is a python script based on [Dash](https://plotly.com/dash/) and [Quickboard](https://github.com/broadinstitute/quickboard/)
to create an interactive dashboard using the outputs of the above WDL. Currently, the board assumes file outputs will exist
for each of the three main components of the workflow. To use, follow the instructions below.

### Dependencies

Follow these steps to ensure your environment is prepared.

1. Make sure you have the Python package `quickboard` installed. It's recommended to pin the version to avoid any breaking
changes in the future. The latest `quickboard==0.3.5` should work with this current version.
2. Copy all the files from the [SVisualizer](SVisualizer) directory to a fresh directory on your local system. 
3. Copy all the files from the WDL outputs into a subdirectory called `wdl_outputs`. You can also use the gathering script
[provided](#optional-data-gathering-script) to make this step easier.

### Usage

Once you've followed all the steps from the [Dependencies](#dependencies), using the script is straightforward. You can
either open and run the `.ipynb` file provided in Jupyter, or run the `.py` version from the command line. Both will
start a local Flask server (hosted on port 8050 by default) you can navigate to in order to interact with the dashboard.
It will be automatically populated using the WDL output data in the `wdl_outputs` directory. The notebook version
has some markdown headers to make it easier to navigate the different blocks of code, if you wish to modify it for your
own purposes. It also makes it easier to see what variables/settings can be toggled at the top of the page.


### Optional Data Gathering Script

Because there are a lot of files output by the workflow, you might want to use the provided `gather_terra_data.py` script
to automatically set up your environment. This works using the Firecloud API to grab your files from a Terra submission
and put them in the right location, combining outputs from separate workflows run under the same submission. 
To use it, you need to:

1. Make sure you have the Firecloud API Python package installed, e.g. `pip install firecloud`.
2. Make sure you have `gsutil` installed and configured with the correct permissions relative to your Terra workspace.
3. Run `python gather_outputs.py <TERRA_NAMESPACE> <TERRA_WORKSPACE_NAME> <SUBMISSION_ID>`.

This will create the `wdl_outputs` directory with the right files inside. You can find the command arguments by navigating
to the Job History tab and clicking on your desired run. You can use e.g. the Terra URL to figure out most of the arguments.
For example, the URL should look like:
```
app.terra.bio/#workspaces/<TERRA_NAMESPACE>/<TERRA_WORKSPACE_NAME>/job_history/<SUBMISSION_ID>
```

The script works by downloading the combined files `.tar` directory output from the workflow, extracting it to a local
temp directory, and then combining these outputs across all successful workflow runs under the same submission. Note
separate "Experiment" label inputs can be used to distinguish parts of the output tables corresponding to different runs.
The script naively concatenates the files across the runs, so doing this for a large number of workflows inside the same
submission could result in large tables, which may not be ideal for the dynamic SVisualizer.