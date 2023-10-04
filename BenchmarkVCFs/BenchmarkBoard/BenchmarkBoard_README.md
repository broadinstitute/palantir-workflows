# BenchmarkBoard

BenchmarkBoard is an interactive visualizer meant to complement the [SimpleBenchmark](../SimpleBenchmark.wdl) WDL tool for benchmarking VCFs. It is easy to setup for use with data output from the WDL, and can be used as-is as a comprehensive solution for analyzing benchmarking data produced by the pipeline, or thought of as a template or inspiration for ways to visualize the large amounts of data produced.

## Quickstart

To use the BenchmarkBoard, you must satisfy the following requirements:
1. Have a Python environment with [Quickboard](https://github.com/broadinstitute/quickboard) installed, e.g. `pip install quickboard`. The current version of this app should use `quickboard==0.3.3`. 
2. Have the files output by the SimpleBenchmark tool saved locally with their default names.
3. Have a copy of the `BenchmarkBoard.py` or `BenchmarkBoard.ipynb` file saved in the *same* directory as the files from the WDL.

Once you satisfy those requirements, you can either run the `.py` file or run all the cells in the `.ipynb` file to start up the dashboard server. You can connect to it by going to `localhost:8050` in any web browser, or clicking the output link if running in Jupyter. If this port is already in use for another application, you can do a find and replace in the Python files for `port=8050` and replace with a free port on your system.

## Overview

When opening the dashboard, you will see various tabs on the top of the screen with can be clicked on to see new pages of plots. These are separated into a few different categories based on the different files produced by the WDL. On the side and below many plots, there are plugins, i.e. buttons or sliders, that can manipulate the data being displayed in the plot. In particular, every page has a sidebar plugin which controls which interval file the data is being stratified by, and what type of variant data is being displayed.

When your data is detected to have multiple samples, i.e. concatenated files across multiple runs of the WDL, the behavior of the board will adapt. Bar and line charts will include "error bars" modeled on a naive normal distribution, and values will be replaced with averages over the `Experiment` groups. 

## Custom Usage

Aside from allowing for user-defined experiment groups to manipulate the multi-sample board, there are a few other options users can configure to change the behavior of the dashboards produced, by editing the top of the Python files themselves.

In the multi-sample board, users can find the following variables to be optionally edited:
* `COVARIATE_X` - (default = None) This can be a string, or list of strings, which are column names the user might want to control for when doing comparisons of statistics like Precision, Recall, etc. When turned on (via providing legal values), this feature will introduce an extra plot in the `Simple Summary` tab called `Stat Correlation Plot`. The y-axis will plot whatever statistics is desired by the user (e.g. Precision), whereas the x-axis can be toggled via radio buttons to be any entry in the `COVARIATE_X` list. A popular use case would be to use user-provided `MEAN_COVERAGE` as a value, as this often correlates with benchmarking performance. This would require the user to either manually edit the files after running the WDL, or to provide them via e.g. the Terra table through the `extra_column` and `extra_column_label` arguments.
* `SIMPLE_VARIANTS_ONLY` - (default = True) This is a boolean which controls which types of variants are available as options when plotting on the `Simple Summary` tab. When `True`, this restricts the options to the usual variant types: SNP, HetSNP, HomVarSNP, INDEL, HetINDEL, and HomVarINDEL. Otherwise, a number of other `bcftools` categories like MA (multiallelic sites), Other, etc. are available. This is toggled off by default, since those types of variants are often much more rare, and their statistics can be more distracting than insightful in many use cases.
* `EXPERIMENT_ORDER` - provide a list of `Experiment` values to control the order these are rendered in the legends in plots that group by `Experiment` value.
* `EXPERIMENT_COLOR_MAP` - provide a dictionary of `Experiment` label to hex color (e.g. `#0055bb`) to override the default plotly palette, and for consistent coloring of `Experiment` groups across plots.
* `SHOW_ERROR_BARS` - if true, show error bars on plots when appropriate. Otherwise globally hide them.