# PALANTIR-WORKFLOWS

[![Test Status](https://github.com/broadinstitute/palantir-workflows/actions/workflows/run_tests.yaml/badge.svg?branch=main)](https://github.com/broadinstitute/palantir-workflows/actions/workflows/run_tests.yaml/badge.svg?branch=main) [![Dockerfiles Status](https://github.com/broadinstitute/palantir-workflows/actions/workflows/test_dockerfiles.yaml/badge.svg?branch=main)](https://github.com/broadinstitute/palantir-workflows/actions/workflows/test_dockerfiles.yaml/badge.svg?branch=main)

Utility workflows used by the DSP's Palantir team.  This repository should be used to manage frequently used utility workflows for the team, and facilitate their use on [Terra](https://app.terra.bio/) through [Dockstore](https://dockstore.org/).

**Remember, this is a public repository, so anything you put in this repo is publicly viewable.**


## Directory Index

Each directory below groups a set of related workflows and has its own README with per-workflow inputs, outputs, and usage notes.

### Benchmarking and evaluation
* [BenchmarkVCFs](BenchmarkVCFs/README.md): Benchmarking short germline variants in VCFs against a truth set using `vcfeval`, plus workflows for finding matching samples, comparing benchmark runs, trio analysis, and a `vcfdist`-based alternative.
* [BenchmarkSVs](BenchmarkSVs/README.md): Benchmarking structural variants using `wittyer`, with an accompanying SVisualizer dashboard.
* [BenchmarkPhasing](BenchmarkPhasing/README.md): Adding read-backed phasing information to a VCF (`whatshap`) and benchmarking phasing between experimental conditions.
* [FunctionalEquivalence](FunctionalEquivalence/README.md): Evaluating whether data processed by different pipelines can be used interchangeably without batch effects.
* [LongReadRNABenchmark](LongReadRNABenchmark/README.md): Comparing long-read RNA isoform discovery tools (Bambu, Cupcake, FLAIR, FLAMES, IsoQuant, IsoSeq, StringTie, TALON) against each other and against a reference annotation.

### Pipelines
* [GlimpseImputationPipeline](GlimpseImputationPipeline/README.md): Low-coverage genotype imputation with GLIMPSE2, including reference-panel splitting, batched imputation, merging, and QC.
* [PRS](PRS/README.md): The polygenic risk score pipeline — scoring, population PCA and ancestry adjustment, aggregation, and QC — plus [ScoreBGE](PRS/ScoreBGE/README.md) for blended genome-exome scoring and [Validation](PRS/Validation/README.md) workflows.
* [gCNV](gCNV/README.md): Germline copy number variant calling with GATK gCNV in cohort and case mode, single-sample filtering, and wrappers that merge calls for Fabric.
* [HPVDeepSeek](HPVDeepSeek/README.md): The HPV DeepSeek duplex-UMI assay — consensus alignment, HPV genotyping, somatic variant calling, integration breakpoints and sublineage assignment, and depth normalization.

### Single-purpose workflows
* [Fingerprints](Fingerprints/README.md): Collecting Picard fingerprint metrics for a sample.
* [HaplotypeMap](HaplotypeMap/README.md): Building the haplotype map files that the fingerprinting workflows consume.
* [JointCalling](JointCalling/README.md): Creating the sample map file used as input to joint genotyping.

### Utilities
* [Utilities](Utilities/README.md): General purpose tooling — [miscellaneous small WDLs](Utilities/WDLs/README.md), [docker images](Utilities/Dockers/README.md), [interval files](Utilities/IntervalFiles/README.md), and [notebooks](Utilities/Notebooks/README.md).
* [MultiQC_Terra](MultiQC_Terra/MultiQC.md): A notebook that gathers QC files from a Terra workspace's submissions and runs MultiQC over them.
* [Scripts/monitoring](Scripts/monitoring/README.md): A script for collecting memory, CPU, and disk usage from a running task.


## Testing Workflows

Automated WDL testing is implemented using [watt](https://github.com/rickymagner/watt).
To add tests, update  [test/watt_config.yml](test/watt_config.yml).
See the watt documentation for usage details.  Tests are run using github actions, controlled by [.github/workflows/run_tests.yaml](.github/workflows/run_tests.yaml).  Automated testing also validates all WDLs found in the repository, regardless of whether tests have been added for the particular WDL, using `womtool validate`.  Automated testing is performed on all PR's, as well as any pushes to the `main` branch.

During testing, WDLs are run on GCP, with the execution bucket `gs://palantir-workflows-test-execution`.  Generally, input and expected output files are stored in `gs://palantir-workflows-test-data`.  The stdout from watt will be visible in the github actions UI, and cromwell log files can be found in a zipped artifact named `cromwell_logs`, also in the github actions UI.   

## Using the Dockstore Github App to Automatically Update Workflows in Dockstore/Terra
Workflows registered in Dockstore can be automatically synced when changes are pushed to this repo by adding their information to `.dockstore.yml`. 
In this way, a change pushed to a branch in this repo can be automatically propagated into any Terra workspaces using the workflow. 
Details can be found at Dockstore's Github App [documentation](https://docs.dockstore.org/en/develop/getting-started/github-apps/github-apps.html).
