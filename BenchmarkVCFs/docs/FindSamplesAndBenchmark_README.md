# FindSamplesAndBenchmark Guide

This guide describes how to use the [FindSamplesAndBenchmark](../FindSamplesAndBenchmark.wdl) WDL. The goal of the workflow is to automatically match samples using fingerprints and then run [BenchmarkVCFs](../BenchmarkVCFs.wdl) on the matching pairs. Usage details can be found below.

## Acknowledgements

This tool was inspired by the original tool with the same name (first written by Megan Shand), and developed subsequently by various hydro.gen members.

## Workflow Overview

The outline is simple and follows three basic steps:
1. Run [MatchFingerprints](/Utilities/WDLs/MatchFingerprints.wdl) across all possible pairs of `base_vcfs` and `query_vcfs`. This task returns just those which are a fingerprint match as pairs in an array.
2. Scatter over the matching pairs and run `BenchmarkVCFs` with the provided parameters.
3. Combine all of the benchmarking outputs into just one table. Individual benchmarking comparisons can be differentiated by the combination of base/query sample names from the VCFs.


## Workflow Inputs

The inputs are mostly the same as those for `BenchmarkVCFs`, so consult the corresponding [documentation](BenchmarkVCFs_README.md) for details. A few fingerprinting inputs are exposed here, which can be explained in the `MatchFingerprints` [documentation](/Utilities/WDLs/MatchFingerprints.wdl). The biggest difference with this workflow is now the `base_vcfs` and `query_vcfs` inputs accept an array of VCFs that will be matched and benchmarked against each other. 



