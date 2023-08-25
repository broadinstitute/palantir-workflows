# BenchmarkBoard

BenchmarkBoard is an interactive visualizer meant to complement the [SimpleBenchmark](../SimpleBenchmark.wdl) WDL tool for 
benchmarking VCFs. It is easy to setup for use with data output from the WDL, and can be thought of as a template or inspiration 
for ways to visualize the large amounts of data produced.

## Quickstart

To use the BenchmarkBoard, you must satisfy the following requirements:
1. Have a Python environment with [Quickboard](https://github.com/broadinstitute/quickboard) installed, e.g. `pip install quickboard`.
2. Have the files output by the SimpleBenchmark tool saved locally.
3. Have a copy of the benchmarkboard-\*.py or benchmarkboard-\*.ipynb file saved in the *same* directory as the files from 
the WDL. These come in two flavors: SS (single-sample mode) or MS (multi-sample mode). If you have just one sample pair to 
analyze, use the file with `-SS`. Otherwise if you have multiple sample pairs and concatenated the outputs from multiple WDL 
runs, use the file with `-MS`.

Once you satisfy those requirements, you can either run the `.py` file or run all the cells in the `.ipynb` file to start 
up the dashboard server. You can connect to it by going to `localhost:8050` in any web browser, or clicking the output 
link if running in Jupyter. If this port is already in use for another application, you can do a find and replace in the 
Python files for `port=8050` and replace with a free port on your system.

## Overview

When opening the dashboard, you will see various tabs on the top of the screen with can be clicked on to see new pages of 
plots. These are separated into a few different categories based on the different files produced by the WDL. On the side 
and below many plots, there are plugins, i.e. buttons or sliders, that can manipulate the data being displayed in the 
plot. In particular, every page has a sidebar plugin which controls which interval file the data is being stratified 
by, and what type of variant data is being displayed.

Although both boards should support multi-sample data tables, the style in which they display data is geared towards the 
two different applications. The single-sample board contains an interactive table on the `Summary Stats` page, which lets 
users query the output summary statistics. The `Callset Sample` sidebar plugin ensures users only see data for one 
sample at a time, with the option to toggle between them with radio buttons (in the case the underlying data tables have 
been modified via concating outputs from multiple WDL runs).

The multi-sample board, on the other hand, is geared towards visualizing groups of VCFs, with the prototypical application 
being a comparison between different experiment groups, coming from the `Experiment` column defined by users either 
in the WDL run, or in a pre-processing step the user is responsible for. Plots here will be colored by this column, 
if present. The same SNP and INDEL tab plots are in this board, but instead represent averages done over the 
`Experiment` groups. Experiment groups to include in visualizations can be toggled using the sidebar checklist, 
instead of having to manually exclude them in each Plotly plot by clicking on the color group each time a plot is regenerated.

## Custom Usage

Aside from allowing for user-defined experiment groups to manipulate the multi-sample board, there are a few other options 
users can configure to change the behavior of the dashboards produced, by editing the top of the Python files themselves.

In the multi-sample board, users can find the following variables to be optionally edited:
* `stat_correlator` - (default = None) This can be a string, or list of strings, which are column names the user might 
want to control for when doing comparisons of statistics like Precision, Recall, etc. When turned on (via providing legal values), 
this feature will introduce an extra plot in the `Simple Summary` tab called `Stat Correlation Plot`. The y-axis will plot whatever 
statistics is desired by the user (e.g. Precision), whereas the x-axis can be toggled via radio buttons to be any entry 
in the `stat_correlator` list. A popular use case would be to use user-provided `MEAN_COVERAGE` as a value, as this 
often correlates with benchmarking performance. This would require the user to either manually edit the files after running 
the WDL, or to provide them via e.g. the Terra table through the `extra_column` and `extra_column_label` arguments.
* `simple_variants_only` - (default = True) This is a boolean which controls which types of variants are available as 
options when plotting on the `Simple Summary` tab. When True, this restricts the options to the usual variant types: 
SNP, HetSNP, HomVarSNP, INDEL, HetINDEL, and HomVarINDEL. Otherwise, a number of other `bcftools` categories like MA 
(multiallelic sites), Other, etc. are available. This is toggled off by default, since those types of variants are often 
much more rare, and their statistics can be more distracting than insightful.