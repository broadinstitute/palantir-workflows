# Docker Files

This directory contains `dockerimage` files and descriptions of what each is intended for. Here is a description of how
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
* Location: `us.gcr.io/broad-dsde-methods/bcftools:v1.0`
* Used By: [SimpleBenchmark](../../BenchmarkVCFs/SimpleBenchmark.wdl)
* Usage: `bcftools [COMMAND]`
* Version Notes:
  * 1.0: Versions are `bcftools` 1.16

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
    
## rtg

* Directory: RTG
* Description: A docker image containing an installation of [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools). 
In particular, this includes `vcfeval`, and all their other VCF manipulation tools. Also comes with Python (pandas) for
convenient data manipulation tools.
* Location: `us.gcr.io/broad-dsde-methods/rtg:v1.0`
* Used By: [SimpleBenchmark](../../BenchmarkVCFs/SimpleBenchmark.wdl)
* Usage: `rtg [COMMAND]`
* Version Notes: 1.0: Versions are `rtg` 3.12.1