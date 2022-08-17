# Docker Files

This directory contains `dockerimage` files and descriptions of what each is intended for. 

## python-data-slim

* Directory: Python-Data-Slim
* Description: A small docker image with essential data processing and Terra packages in python. Comes with `python3`, along 
with the libraries `pandas`, `numpy`, `scipy`, `firecloud`, `fsspec` (a `firecloud` optional dependency), and `gcsfs` 
(for accessing Google bucket files) explicitly installed.
* Location: `us.gcr.io/broad-dsde-methods/python-data-slim`
* Used By: [FindSamplesAndBenchmark](../../BenchmarkVCFs/FindSamplesAndBenchmark.wdl), 
  [CollectBenchmarkSucceeded](../WDLs/CollectBenchmarkSucceeded.wdl)
* Version Notes: 
    * 1.0: Versions are `python3` 3.9.9, `pandas` 1.3.4, `numpy` 1.21.4, `scipy` 1.7.2, `firecloud` 0.16.32, 
    `fsspec` 2022.7.1, `gcsfs` 2022.7.1.
