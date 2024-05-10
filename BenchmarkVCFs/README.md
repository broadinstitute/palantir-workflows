# Benchmarking VCFs

This directory houses the hydro.gen team's collection of workflows used in benchmarking short germline variants in VCF format. At the center is `BenchmarkVCFs`, which is a dependency for multiple other workflows. Full details on how to use these and interpret their results can be found below. We provide a quick overview, with links to the full documentation.

* [BenchmarkVCFs](docs/BenchmarkVCFs_README.md): Used to generate classifier stats (precision, recall, etc) between a baseline "truth" VCF and a query VCF to measure against it. This uses the tool `vcfeval` developed by [RTG](https://github.com/RealTimeGenomics/rtg-tools) to accurately compare variants according to their generated haplotypes from reference to control for representation mismatches. After acquiring "TP", etc. labels in this way, stats are collected over various subsets controlled by region files, types of variants, and more using `bcftools`. 
* [FindSamplesAndBenchmark](docs/FindSamplesAndBenchmark_README.md): Given two collections of VCFs, fingerprinting is run to match samples across the inputs and run `BenchmarkVCFs` on each matching pair. This can be useful for comparing a cohort with known truth samples to easily pick out the corresponding samples without worrying about name matches.
* [BenchmarkAndCompareVCFs](docs/CompareBenchmarks_README.md): Runs benchmarking across samples, and then creates a table of differences between statistics collected for the different samples. These can be later uploaded to a Google sheet using the provided `ExportToGoogleSheets` notebook.

We also have an experimental workflow [VcfdistBenchmark](docs/VcfdistBenchmark_README.md) written by a previous summer intern to test using the `vcfdist` engine instead of `vcfeval`.
