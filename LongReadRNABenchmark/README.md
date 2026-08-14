# Long Read RNA Isoform Discovery Benchmark

Long read RNA sequencing (PacBio, ONT) produces reads that are often long enough to span a full transcript, so an aligned read can be interpreted directly as a candidate isoform. *Isoform discovery* (also called transcript model reconstruction) is the step that turns a BAM of aligned long reads into a set of transcript models, i.e. a GTF/GFF. Many tools do this, some of them guided by an existing reference annotation and some of them entirely *de novo*, and they disagree with each other quite a bit. The WDLs in this directory run a collection of those tools over the same aligned BAM and score their output against a known truth set, so that they can be compared on equal footing.

The top level workflow, [IsoformDiscoveryBenchmark](#isoformdiscoverybenchmark), fans a single input BAM out to eight per-tool sub-workflows (IsoQuant, StringTie, Bambu, Flair, TALON, IsoSeq, FLAMES, and Cupcake; IsoQuant and StringTie are each run twice, once with and once without a reference annotation), collects the GTF/GFF each one produces, and then evaluates them with [gffcompare](https://github.com/gpertea/gffcompare). Three separate comparisons are made:

* **Reduced annotation analysis.** The tools that take an annotation are given a *reduced* annotation, from which some of the transcripts that are actually expressed in the dataset have been removed. Each tool's GTF is compared against the full set of expressed transcripts (`expressedGTF`) and against the retained subset (`expressedKeptGTF`). Transcripts that are expressed but were removed from the reduced annotation are the "novel" truth set: recovering them is what the tool is being scored on. Sensitivity, precision, and F1 for novel transcripts, plus raw TP/FP/FN counts for novel and known transcripts, are reported per tool.
* **De novo (tool concordance) analysis.** All of the reduced-annotation GTFs are compared against each other in a single gffcompare run, and each transcript model is bucketed by how many tools found it: found by all tools (`Reliable`), all but one (`Almost_Reliable`), several (`Mult_Pred`), only one (`Unique`), or missed by this tool alone (`Missed`). Reported separately for known and novel transcripts.
* **Reference free analysis.** The GTFs produced without an annotation are each compared directly against `expressedGTF` with `gffcompare -r`, and the transcript-level sensitivity/precision reported by gffcompare (plus a derived F1) is tabulated per tool.

Each of the three summaries is written as a TSV and rendered as a bar-chart PNG.

The per-tool sub-workflows are also registered individually and can be run on their own if all you want is one tool's transcript models.

Contents:
- [IsoformDiscoveryBenchmark](#isoformdiscoverybenchmark): run all isoform discovery tools on one BAM and compare their transcript models
- [Per-tool workflows](#per-tool-workflows): [Bambu](#bambu), [Cupcake](#cupcake), [Flair](#flair), [Flames](#flames), [IsoQuant](#isoquant), [IsoSeq](#isoseq), [StringTie](#stringtie), [Talon](#talon)
- [Docker images](#docker-images)
- [Dockstore registration and tests](#dockstore-registration-and-tests)

Note on naming: the file is [`IsoformDiscoveryBenchmark.wdl`](IsoformDiscoveryBenchmark.wdl) and the Dockstore entry is `IsoformDiscoveryBenchmark`, but the WDL `workflow` declaration inside it is named `LongReadRNABenchmark`. That is the name you will see in Cromwell/Terra job history.

## IsoformDiscoveryBenchmark

### Summary

The workflow takes a single BAM of long RNA reads already aligned to the reference genome (it does not do the alignment itself; some tools re-derive FASTQ from the BAM with `samtools fastq` and realign internally), plus a reference genome, a reference annotation, and the two "expressed" truth GTFs. It then runs, in parallel:

1. **[IsoQuant](#isoquant)** with `referenceAnnotation`, and **IsoQuantReferenceFree** with no annotation.
2. **[StringTie](#stringtie)** with `referenceAnnotation`, and **StringTieReferenceFree** with no annotation.
3. **[Bambu](#bambu)**, **[Flair](#flair)**, **[Talon](#talon)**, and **[Flames](#flames)**, all annotation-guided.
4. **[IsoSeq](#isoseq)** and **[Cupcake](#cupcake)**, both annotation-free.

The resulting GTF/GFFs are then collected into two lists, which the evaluation steps use:

* `gtfListReduced` = IsoQuant, StringTie, Bambu, Flair, Talon, Flames (the annotation-guided runs), labelled by `toolNamesReduced`.
* `gtfListReferenceFree` = IsoQuantReferenceFree, StringTieReferenceFree, IsoSeq, Cupcake, labelled by `toolNamesReferenceFree`.

**These lists and their name arrays are hardcoded in the workflow body and must stay in the same order.** There is an explicit comment to this effect in the source: if they get out of sync you will generally not get an error, you will get silently wrong results. The `numTools = 6` default on the `SummarizeDenovoAnalysis` task likewise matches the length of `gtfListReduced`.

Evaluation then proceeds as:

5. `GffCompareTrack` (scattered over the six reduced-annotation GTFs): runs `gffcompare -o {datasetName}.{toolName} {expressedGTF} {expressedKeptGTF} {toolGTF}` and keeps the `.tracking` file. Because the truth GTFs occupy the first two query columns, the tracking file tells you for every transcript whether it is expressed, whether it survived into the reduced annotation, and whether the tool found it.
6. `GffCompareTrackDenovo`: one `gffcompare -o {datasetName}.denovo {expressedKeptGTF} {toolGTFs...}` run across all six tool GTFs at once, for the tool-concordance analysis.
7. `ReferenceFreeAnalysis` (scattered over the four annotation-free GTFs): `gffcompare -r {expressedGTF} -o {basename}.reffree {inputGTF}`, keeping the gffcompare stats output.
8. `SummarizeAnalysis`, `SummarizeDenovoAnalysis`, and `SummarizeReferenceFreeAnalysis`: python scripts (baked into the custom docker image, see [Docker images](#docker-images)) that parse the gffcompare output into TSV summary tables.
9. `PlotAnalysisSummary` (run twice, for the reduced and reference-free summaries) and `PlotDenovoAnalysisSummary` (run twice, for known and novel): matplotlib/seaborn bar charts of the summary tables.

Every per-tool task also writes a `monitoring.log` (resource usage over time, from a monitoring script staged out of a Google bucket). These logs are outputs of the per-tool workflows but are *not* surfaced as outputs of the top level benchmark, so to see them you need to dig into the call outputs in the job history.

### Inputs

* `File inputBAM`: long RNA reads aligned to `referenceGenome`.
* `File inputBAMIndex`: index for `inputBAM`.
* `File referenceGenome`: reference genome FASTA.
* `File referenceGenomeIndex`: index for `referenceGenome`.
* `File referenceAnnotation`: the annotation handed to the annotation-guided tools (IsoQuant, StringTie, Bambu, Flair, Talon, Flames). For the intended benchmark this is a *reduced* annotation, i.e. one with some genuinely expressed transcripts removed, so that the tools have something to rediscover.
* `File expressedGTF`: truth set of transcripts actually expressed in this dataset. Used as the first query to `GffCompareTrack` and as the `-r` reference for the reference-free analysis.
* `File expressedKeptGTF`: the subset of `expressedGTF` retained in the reduced annotation. Transcripts in `expressedGTF` but not in `expressedKeptGTF` are treated as novel for scoring purposes.
* `String datasetName`: label for this dataset. Used in every output filename, as the sample label passed to IsoQuant, and as the TALON build/annotation/dataset name.
* `String dataType`: sequencing data type. Passed straight through to IsoQuant's `--data_type` argument and written into the TALON config CSV as the platform field, so it must be a value IsoQuant accepts (see the [IsoQuant documentation](https://github.com/ablab/IsoQuant)).

Per-task runtime attributes (`cpu`, `numThreads`, `memoryGB`, `diskSizeGB`, `docker`, `monitoringScript`) are not exposed at the top level; they are defaulted on the individual tasks and can only be overridden by editing the WDL or, in Terra/Cromwell, by supplying fully-qualified call inputs.

### Outputs

* `File analysisSummary`: `{datasetName}_analysis_summary.tsv`, one row per reduced-annotation tool with columns `Tool`, `Sensitivity(Novel)`, `Precision(Novel)`, `F1-Score(Novel)`, `TP_Novel`, `FP_Novel`, `FN_Novel`, `TP_Known`, `FN_Known`.
* `File analysisSummaryReferenceFree`: `{datasetName}_analysis_summary_reffree.tsv`, one row per annotation-free tool with columns `Tool`, `Sensitivity`, `Precision`, `F1-Score` (transcript-level, taken from gffcompare and rescaled to 0-1).
* `File analysisSummaryDenovoKnown`: `{datasetName}_denovo_analysis_summary_known.tsv`, a table with one column per tool and rows `Tool_Names`, `Reliable_Transcripts`, `Almost_Reliable_Transcripts`, `Mult_Pred`, `Unique`, `Missed`, restricted to transcripts present in `expressedKeptGTF`. `Missed` is written as a negative number so that it plots below the axis.
* `File analysisSummaryDenovoNovel`: same table as above, for transcripts *not* present in `expressedKeptGTF`.
* `File analysisSummaryPlot`: `{datasetName}_analysis_summary_reduced.png`, three stacked bar panels (sensitivity, precision, F1) over the reduced-annotation tools.
* `File referenceFreeAnalysisSummaryPlot`: `{datasetName}_analysis_summary_reffree.png`, same three panels over the annotation-free tools.
* `File denovoAnalysisSummaryPlotKnown`: `{datasetName}_analysis_summary_denovo_known.png`.
* `File denovoAnalysisSummaryPlotNovel`: `{datasetName}_analysis_summary_denovo_novel.png`.

## Per-tool workflows

Each of these is a standalone, individually registered workflow wrapping a single task. They all share the same shape: one long-running task, a `monitoring.log` output, and hardcoded resource defaults (`cpu`, `numThreads`, `memoryGB`, `diskSizeGB`) and a pinned `docker` digest on the task. Each also takes a `monitoringScript` input, defaulting to `gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh`, which is run in the background to record resource usage.

Unless noted otherwise, `inputBAM` is expected to be long RNA reads already aligned to `referenceGenome`.

### Bambu

[`Bambu.wdl`](Bambu.wdl). Runs [bambu](https://bioconductor.org/packages/bambu) (an R/Bioconductor package) via an inline `Rscript`: `prepareAnnotations()` on the reference annotation, then `bambu()` on the BAM, then `writeBambuOutput()`. The task then post-processes bambu's `extended_annotations.gtf` down to only those transcripts with a read count of at least 1 in `counts_transcript.txt`.

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `referenceAnnotation` (required), `datasetName`.
* Outputs: `bambuGTF` (`Bambu_out/Bambu_out_{datasetName}.gtf`, the expression-filtered extended annotation), `bambuCounts` (the filtered per-transcript count table), `monitoringLog`.

### Cupcake

[`Cupcake.wdl`](Cupcake.wdl). Converts the BAM back to FASTQ with `samtools fastq`, removes duplicate FASTQ records with the bundled `remove_fastq_duplicates.py` helper, then runs [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)'s `collapse_isoforms_by_sam.py` to collapse the aligned reads into isoforms. Annotation-free. Note the collapse step is run with `--cpus 1` regardless of the requested CPU count.

* Inputs: `inputBAM`, `inputBAMIndex`, `datasetName`.
* Outputs: `cupcakeGFF` (`Cupcake_out_{datasetName}.collapsed.gff`), `monitoringLog`.

### Flair

[`Flair.wdl`](Flair.wdl). Runs [FLAIR](https://github.com/BrooksLabUCSC/flair) from the Docker Hub image `brookslab/flair`: `samtools fastq` to recover reads, `bam2Bed12` to convert the BAM to BED12, then `flair correct` (splice-site correction against the genome and annotation) followed by `flair collapse` (isoform collapsing).

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `referenceAnnotation` (required), `datasetName`.
* Outputs: `flairGTF` (`Flair_out_{datasetName}.isoforms.gtf`), `monitoringLog`.

### Flames

[`Flames.wdl`](Flames.wdl). Runs [FLAMES](https://github.com/LuyiTian/FLAMES)' `python/bulk_long_pipeline.py`, giving it both the original BAM (`--inbam`) and a FASTQ directory produced from it with `samtools fastq`. Note that FLAMES expects the annotation as GFF3 (`--gff3`).

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `referenceAnnotation` (required). Note this is the only per-tool workflow that does **not** take `datasetName`, so its output filename is fixed.
* Outputs: `flamesGFF` (`isoform_annotated.gff3`), `monitoringLog`.

### IsoQuant

[`IsoQuant.wdl`](IsoQuant.wdl). Runs [IsoQuant](https://github.com/ablab/IsoQuant)'s `isoquant.py` on the BAM. The reference annotation is optional: when supplied, the task passes `--genedb` plus `--complete_genedb` and writes to `IsoQuant_out_{datasetName}`; when omitted, it runs annotation-free and writes to `IsoQuant_denovo_out_{datasetName}`. This is how the top level workflow gets both an annotation-guided and a reference-free IsoQuant run out of one wrapper.

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `referenceAnnotation` (**optional**), `datasetName`, `dataType`.
* Outputs: `isoQuantGTF` (`{outputPrefix}/{datasetName}/{datasetName}.transcript_models.gtf`), `monitoringLog`.

### IsoSeq

[`IsoSeq.wdl`](IsoSeq.wdl). Converts the BAM to FASTQ, realigns with [pbmm2](https://github.com/PacificBiosciences/pbmm2) using the `ISOSEQ` preset, and collapses redundant transcripts with [isoseq3](https://github.com/PacificBiosciences/IsoSeq) `collapse`. Annotation-free. Because it realigns from scratch, this is one of the more expensive tasks in the benchmark.

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `datasetName`.
* Outputs: `isoSeqGFF` (`IsoSeq_out_{datasetName}.gff`), `monitoringLog`.

### StringTie

[`StringTie.wdl`](StringTie.wdl). Runs [StringTie](http://ccb.jhu.edu/software/stringtie/) in long-read mode (`-L`). As with IsoQuant, the annotation is optional: with it, StringTie runs in reference-guided mode (`-G`) and the output is named `StringTie_out_{datasetName}.gtf`; without it, the output is named `StringTie_denovo_out_{datasetName}.gtf`.

* Inputs: `inputBAM`, `referenceAnnotation` (**optional**), `datasetName`. Note this workflow does not take a BAM index or the reference genome.
* Outputs: `stringTieGTF`, `monitoringLog`.

### Talon

[`Talon.wdl`](Talon.wdl). Runs the full [TALON](https://github.com/mortazavilab/TALON) pipeline: `talon_label_reads` (internal priming labels), `samtools calmd` to add MD tags, `talon_initialize_database` from the reference annotation, a generated single-line config CSV, `talon` itself, `talon_filter_transcripts`, and finally `talon_create_GTF` against the filtered whitelist. `datasetName` is reused as the TALON build name, annotation name, and dataset name; `dataType` is the platform field of the config CSV.

* Inputs: `inputBAM`, `inputBAMIndex`, `referenceGenome`, `referenceGenomeIndex`, `referenceAnnotation` (required), `datasetName`, `dataType`.
* Outputs: `talonGTF` (`Talon_out_{datasetName}_talon.gtf`), `monitoringLog`.

## Docker images

Every task pins its image by digest (`@sha256:...`). With the exception of Flair, which uses the public Docker Hub image `brookslab/flair`, all images live under `us.gcr.io/broad-dsde-methods/kockan/` and are built from the Dockerfiles in this directory.

**There are no build/push scripts checked in for these images.** To rebuild one, build the directory and push it, then update the pinned digest in the corresponding WDL, e.g.:

```
cd isoquant_docker
docker build -t us.gcr.io/broad-dsde-methods/kockan/isoquant:<tag> .
docker push us.gcr.io/broad-dsde-methods/kockan/isoquant:<tag>
```

Note that not all of these builds are reproducible: the FLAMES and cDNA_Cupcake images `git clone` the tool's default branch rather than a tagged release, and the IsoSeq image installs from bioconda without version pins, so rebuilding those will not necessarily reproduce the currently pinned image. The rest install from pinned release tarballs or version-pinned packages.

| Directory | Used by | Contents |
| --- | --- | --- |
| [`bambu_docker`](bambu_docker) | [Bambu](#bambu) | Ubuntu 20.04 + R (installed from the RStudio CDN .deb) + BiocManager + the `bambu` Bioconductor package. |
| [`cupcake_docker`](cupcake_docker) | [Cupcake](#cupcake) | Ubuntu 20.04, samtools, BioPython, bx-python, Cython, a deliberately pinned old numpy (cDNA_Cupcake fails on newer numpy), cDNA_Cupcake cloned from GitHub and installed, plus the local helper script `remove_fastq_duplicates.py` copied to `/usr/local/src`. Sets `SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL` to work around cupcake depending on the deprecated `sklearn` shim package. |
| [`flames_docker`](flames_docker) | [Flames](#flames) | Ubuntu 20.04, minimap2, samtools, pysam/numpy/editdistance, and FLAMES cloned to `/usr/local/src/FLAMES`. |
| [`gffcompare_docker`](gffcompare_docker) | all evaluation tasks in [`IsoformDiscoveryBenchmarkTasks.wdl`](IsoformDiscoveryBenchmarkTasks.wdl) | Ubuntu 20.04 with gffcompare built from a release tarball and installed to `/usr/local/bin`. |
| [`isoquant_docker`](isoquant_docker) | [IsoQuant](#isoquant) | Ubuntu 20.04, minimap2, samtools, and IsoQuant unpacked to `/usr/local/src/IsoQuant-<version>` with its requirements installed. The WDL invokes the versioned path directly, so bumping the IsoQuant version here requires editing the path in `IsoQuant.wdl` too. |
| [`isoseq_docker`](isoseq_docker) | [IsoSeq](#isoseq) | `condaforge/mambaforge` base with samtools, pbmm2, and isoseq3 installed from bioconda. |
| [`lr_isoform_custom_docker`](lr_isoform_custom_docker) | the summarize/plot tasks in [`IsoformDiscoveryBenchmarkTasks.wdl`](IsoformDiscoveryBenchmarkTasks.wdl) | Ubuntu 20.04 with matplotlib, numpy, pandas, and seaborn, plus the five analysis scripts copied to `/usr/local/src`: `summarize_analysis.py`, `summarize_reffree_analysis.py`, `summarize_denovo_analysis.py`, `plot_analysis_summary.py`, `plot_denovo_analysis_summary.py`. If you change any of those scripts you must rebuild and repin this image, since the WDL runs them out of the image rather than localizing them. |
| [`stringtie_docker`](stringtie_docker) | [StringTie](#stringtie) | Ubuntu 20.04, minimap2, samtools, and StringTie built from source. |
| [`talon_docker`](talon_docker) | [Talon](#talon) | Ubuntu 20.04, samtools, bedtools (static binary), cython/pytest, and TALON installed from a release tarball. |

There is no Dockerfile here for Flair; that workflow uses the upstream `brookslab/flair` image directly.

## Dockstore registration and tests

All nine workflows in this directory are registered in the repo-level [`.dockstore.yml`](../.dockstore.yml): `IsoformDiscoveryBenchmark` (pointing at `IsoformDiscoveryBenchmark.wdl`), plus `IsoQuant`, `StringTie`, `Bambu`, `Flair`, `Talon`, `IsoSeq`, `Flames`, and `Cupcake`. `IsoformDiscoveryBenchmarkTasks.wdl` is a task library imported by the top level workflow and is not registered separately.

None of these workflows are currently covered by the repo's automated WATT tests: there is no entry for them in [`test/watt_config.yml`](../test/watt_config.yml) and no input JSONs under `test/`. Changes here should be validated by running the workflow manually.
