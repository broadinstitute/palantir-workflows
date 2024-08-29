# CompareBenchmarks

The purpose of this workflow is to create a table to facilitate a comparison of precision and sensitivity between different configurations (e.g. pipeline versions, chemistry changes) of samples vs. truth that have been obtained with the `BenchmarkVCFs` workflow.

After the `CompareBenchmarks` workflow is run, you can export the generated CSV table into a Google Sheets spreadsheet with automatic formatting using the `ExportToGoogleSheets` Colab notebook described below.

## Workflow

### Input

The input is provided using 3 arrays: `sample_ids`, `configurations` and `benchmark_summaries`. For each sample/configuration combination one entry in these arrays is required. For example, if NA12878 called with GATK3 should be compared with NA12878 called with GATK4, the input should look like this:

```
sample_ids          = ["NA12878", "NA12878"]
configurations      = ["GATK3", "GATK4"]
benchmark_summaries = [{path_to_NA12878_GATK3_vs_truth}, {path_to_NA12878_GATK4_vs_truth}]
```

It is possible to provide more than two different configurations and more than one sample, but **the order between the corresponding entries in the arrays must be consistent**. That means that the lengths of these arrays will always be equal. An example for comparing three configurations to each other using two samples is provided in the test files (see `/test/CompareBenchmarks/`).

- `Array[String] sample_ids`: The names of one or more different samples. The comparisons will be made between configurations for each sample individually.
- `Array[String] configurations`: The labels for the different configurations that should be compared to each other.
- `Array[File] benchmark_summaries`: The output `SimpleSummary.tsv` files from the `BenchmarkVCFs` workflow.
- `Array[String]? stratifiers`: This input requires the same labels as the `stratifier_labels` input that has been passed to `BenchmarkVCFs`.
- `Boolean include_counts = false`: If set to false, the resulting metrics will be Sensitivity, Precision and F-Measure. If set to true, the output table will also include the number of TP, FP and FN variants for each stratifier.
- `Array[String]? order_of_samples`: If multiple different sample names are provided you can specify the order of those samples in the resulting table. Here, each sample name only has to be specified once, not once for each input VCF. If not specified, the order will be determined by the input files.
- `Array[String]? order_of_configurations`: This input determines the order of the configurations in the resulting table. Just as above, each configuration only has to be specified once, not once for each input VCF. If not specified, the order will be determined by the input files.
- `Array[Int]? deltas`: If specified, the output table will contain columns that compares the results of different configurations to each other. Each delta column is defined by two entries in the array (referenced by zero-based index). If, for example there are three configurations A, B and C and you want to compare configurations B to A and C to A, provide the following data: `[1, 0, 2, 0]`. Note that you will want to define the `order_of_configurations` in this case to make sure that the indices refer to the correct configurations.
- `Int? mem_gb`: Optional input overriding the default memory.
- `Int? preemptible`: Optional input overriding the default number of preemptible attempts.


## Notebook for Exporting Results to Google Sheets

This directory also includes a Google Colab notebook `ExportToGoogleSheets.ipynb`. In order to authenticate and access Google Sheets, **this notebook must be run in Google Colab** and cannot be run in a Jupyter environment. The very first cell of the notebook contains a button to open it in Google Colab.

This notebook accepts a path to a comparison CSV file that has been created with this workflow and adds it as a new sheet (tab) to an existing spreadsheet and automatically formats that sheet including borders and conditional coloring of the delta columns.