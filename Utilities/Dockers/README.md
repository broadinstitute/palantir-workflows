# Docker Files

This directory contains `Dockerfile` files and descriptions of what each is intended for. Here is a description of how
the fields are populated:
* Directory: where to find the Dockerfile, under `palantir-workflows/Utilities/Dockers/`.
* Description: description of the packages installed on the docker.
* Location: where a pre-built copy of the docker is stored.
* Used By: links to all WDLs in this repo which use this docker.
* Usage: command line notation for accessing primary installed packages within the docker.
* Version Notes: notes on which versions of major packages are installed on the image.

## bcftools

* Directory: BCFTools
* Description: A docker image containing an installation of [bcftools](https://samtools.github.io/bcftools/bcftools.html).
Also includes some minimal Python tools (pandas) for data processing.
* Location: `us.gcr.io/broad-dsde-methods/bcftools:v1.1`
* Used By: [SimpleBenchmark](../../BenchmarkVCFs/SimpleBenchmark.wdl)
* Usage: `bcftools [COMMAND]`
* Version Notes:
  * 1.0: Versions are `bcftools` 1.16.
  * 1.1: Added the `python-is-python3` package so that `python` is a valid command.

## bedtools

* Directory: Bedtools
* Description: A docker image containing an installation of [bedtools](https://bedtools.readthedocs.io/en/latest/).
* Location: `us.gcr.io/broad-dsde-methods/bedtools:v1.0`
* Used By: [BenchmarkPhasing](../../BenchmarkPhasing/BenchmarkPhasing.wdl)
* Usage: `bedtools [COMMAND]`
* Version Notes:
  * 1.0: Versions are `bedtools` 2.30.0

## samtools

* Directory: Samtools
* Description: A docker image containing an installation of [samtools](https://github.com/samtools/samtools). Also
includes some minimal Python tools (pandas) for data processing.
* Location: `us.gcr.io/broad-dsde-methods/samtools:v1`
* Used By: [ComputeIntervalBamStats](../IntervalFiles/ComputeIntervalBamStats.wdl)
* Usage: `samtools [COMMAND]`
* Version Notes:
  * 1.0: Versions are `samtools` 1.16.1.

## pysam

* Directory: Pysam
* Description: A docker image containing Python with the [Pysam](https://pysam.readthedocs.io/en/latest/api.html) package.
* Location: `us.gcr.io/broad-dsde-methods/pysam:v1`
* Used By: ...
* Usage: `python <<CODE [CODE HERE] CODE`
* Version Notes:
  * 1.0: Versions are `python3` 3.11.2, `pysam` 0.20.0.

## python-data-slim

* Directory: Python-Data-Slim
* Description: A small docker image with essential data processing and Terra packages in python. Comes with `python3`, along 
with the libraries `pandas`, `numpy`, `scipy`, `firecloud`, `fsspec` (a `firecloud` optional dependency), and `gcsfs` 
(for accessing Google bucket files) explicitly installed.
* Location: `us.gcr.io/broad-dsde-methods/python-data-slim`
* Used By: [FindSamplesAndBenchmark](../../BenchmarkVCFs/FindSamplesAndBenchmark.wdl), 
  [CollectBenchmarkSucceeded](../WDLs/CollectBenchmarkSucceeded.wdl)
* Usage: `python <<CODE [CODE HERE] CODE`
* Version Notes: 
  * 1.0: Versions are `python3` 3.9.9, `pandas` 1.3.4, `numpy` 1.21.4, `scipy` 1.7.2, `firecloud` 0.16.32, 
  `fsspec` 2022.7.1, `gcsfs` 2022.7.1.

## python-data-slim-plots

* Directory: Python-Data-Slim-Plots
* Description: Based on `python-data-slim`, this image supports plotting with [matplotlib](https://matplotlib.org/) and [Plotly](https://plotly.com/)
* Location: `us.gcr.io/broad-dsde-methods/python-data-slim-plots`
* Used By: [FunctionalEquivalence](../../FunctionalEquivalence/FunctionalEquivalence.wdl)
* Version Notes:
  * 1.0: Versions are `matplotlib` 3.5.3, `plotly` 5.10.0.
    
## rtg

* Directory: RTG
* Description: A docker image containing an installation of [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools). 
In particular, this includes `vcfeval`, and all their other VCF manipulation tools. Also comes with Python (pandas) for
convenient data manipulation tools.
* Location: `us.gcr.io/broad-dsde-methods/rtg:v1.0`
* Used By: [SimpleBenchmark](../../BenchmarkVCFs/SimpleBenchmark.wdl)
* Usage: `rtg [COMMAND]`
* Version Notes: 
  * 1.0: Versions are `rtg` 3.12.1

## whatshap

* Directory: WhatsHap
* Description: A docker image containing an installation of [WhatsHap](https://whatshap.readthedocs.io/en/latest/).
* Location: `us.gcr.io/broad-dsde-methods/whatshap:v1`
* Used By: [BenchmarkPhasing](../../BenchmarkPhasing/BenchmarkPhasing.wdl), [PhaseVCF](../../BenchmarkPhasing/PhaseVCF.wdl)
* Usage: `whatshap [COMMAND]`
* Version Notes:
  * 1.0: Versions are `whatshap` 1.7


## vcfdist

* Directory: Vcfdist
* Description: A docker image containing an installation of [vcfdist](https://github.com/TimD1/vcfdist).
* Location: ``
* Used By: ``
* Usage: ``
* Version Notes: