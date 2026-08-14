# Utilities

This directory collects general purpose tooling that supports the rest of the repo: small standalone workflows, the
docker images those workflows run on, reference interval files, and one-off analysis notebooks. Each subdirectory has
its own README with the details; this page is just an index.

* [WDLs](WDLs/README.md): a collection of miscellaneous WDLs useful for small tasks, e.g. annotating a VCF, indexing a
  CRAM, converting interval lists to beds, matching fingerprints, or collecting RNA-seq metrics.
* [IntervalFiles](IntervalFiles/README.md): documentation and bucket mirrors for the popular interval files (mostly from
  NIST's GIAB genome stratifications) used to restrict analysis to subsets of the genome, along with overlap statistics
  between them and a workflow for computing BAM statistics over them.
* [Dockers](Dockers/README.md): the `Dockerfile`s used by workflows in this repo, with a description of the packages
  installed on each, where a pre-built copy lives, and which WDLs use it.
* [Notebooks](Notebooks/README.md): Jupyter notebooks for one-off analyses or computational tools.
