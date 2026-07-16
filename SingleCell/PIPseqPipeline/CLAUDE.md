# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow DSL2 pipeline for processing PIPseq single-cell RNA-seq data: it runs DRAGEN scRNA per subsample, then produces per-subsample and supersample-level QC reports, with optional CRISPR guide assignment (CRISPAT). It's part of the larger `palantir-workflows` monorepo (see `../../.github/copilot-instructions.md` for repo-wide conventions) but is developed and run from this directory. The deployment target is Illumina Connected Analytics (ICA) — DRAGEN runs on ICA's FPGA-preset pods (see `pod annotation:` lines in `modules/dragen_scrna.nf`), and `export_pipeline_to_ica.py` (plus `nextflow_schema.json`, which drives the ICA-rendered input form) exist to publish this pipeline into ICA.

**README.md and docker/README.md are stale.** They describe an older CLI surface (`--samplesheet`, `--data_dir`, `-profile local/awsbatch/gcp/cluster/docker`) and multiple per-process Docker images that no longer match the code. There is currently no `profiles {}` block in `nextflow.config` at all, and the pipeline takes `--fastq_list` (not `--samplesheet`). Trust `main.nf`, `nextflow_schema.json`, and `modules/*.nf` over the READMEs when they conflict.

## Terminology

- **Subsample**: an individual sequencing unit (e.g. one well/technical replicate), identified by `RGSM` in the fastq list.
- **Supersample**: a logical grouping of subsamples analyzed together, identified by `--supersample_id` / `--supersample_basename`.

## Running / testing

There's no unit test suite; correctness is verified by running the pipeline (or Nextflow's `-stub-run` mode) end to end.

```bash
# Stub run (no real DRAGEN/CRISPAT execution — every module has a `stub:` block that
# just touches placeholder output files so you can validate pipeline wiring quickly)
cd stub_test
nextflow run ../main.nf -stub-run -params-file stub_inputs/pipeline_input.json

# Real run
nextflow run main.nf --num_input_cells <int> --fastq_list <fastq_list.csv> \
  --supersample_id <id> --supersample_basename <name> \
  --min_valid_guides <int> --max_valid_guides <int> \
  --ref_tar <ref.tar> --annotation_file <annotation.gtf> \
  --dragen_container <image> --qc_container <image> \
  --concatenate_cpus <n> --concatenate_memory_gb <n> --dragen_scratch_tb <n>
```

Required params are enforced by hand-rolled checks at the top of the `workflow` block in `main.nf` (not by Nextflow schema validation at runtime) — if you add a new required param, add a corresponding `if (!params.x) { log.error ...; exit 1 }` check there, and also update `nextflow_schema.json` (used by ICA to render the input form) and the help text in `helpMessage()`.

`stub_test/stub_inputs/` has several ready-made `pipeline_input*.json` variants (with/without feature barcodes, with cell hashing) for exercising different branches of the fastq-list parsing logic.

## Architecture

`main.nf` is a single linear `workflow {}` block (no subworkflows). Data flows through modules in `modules/`:

1. **Parse `--fastq_list`** (CSV with `RGID, RGSM, RGTY, Read1File, Read2File`). Rows are grouped by `RGSM` (subsample). `RGTY` is `'expression'`, `'feature'`, or `'hashing'`; feature/hashing RGIDs are collected per-subsample into comma-joined lists (`scrna_feature_barcode_groups`, `scrna_hto_barcode_groups`) that get passed straight through to DRAGEN's `--scrna-direct-capture-barcode-groups`/`--scrna-feature-barcode-groups`/`--scrna-hto-barcode-groups` args.
2. **`DRAGEN_SCRNA`** (`modules/dragen_scrna.nf`) — one invocation per subsample. Emits a flat `dragen_output/*` file glob (metrics CSV, barcode summary TSV, filtered matrix/barcodes/features). `main.nf` re-derives structure from these by parsing filenames (`file_name.tokenize('.')[0]` for subsample id, suffix matching for file type) and `groupTuple`-ing back into per-subsample tuples — if you change DRAGEN's `--output-file-prefix` or output filenames, this parsing breaks silently.
3. **`GENERATE_REPORT_DATA`** (`modules/generate_report_data.nf`, wraps `bin/generate_report_data.py`) — runs per subsample, always executes regardless of guide assignment.
4. **`CONCATENATE`** (`modules/concatenate.nf`, wraps `bin/concatenate_samples.py`) — always runs (independent of `run_guide_assignment`); merges all subsamples' matrices into one supersample-level `.h5ad` plus a CRISPR-features-only `.crispr.h5ad`. Handles the single-subsample case automatically (no special-casing needed by the caller).
5. **`GUIDE_ASSIGNMENT`** (`modules/guide_assignment.nf`, wraps `bin/run_guide_assignment.py`) — only runs if `params.run_guide_assignment` (default `true`); uses CRISPAT's Poisson-Gaussian mixture model on `CONCATENATE`'s CRISPR AnnData output. When disabled, downstream code substitutes a `Channel.of(file('NO_FILE'))` sentinel — modules check `guide_assignments.name != 'NO_FILE'` rather than checking optionality via Nextflow's `optional` input mechanism.
6. **`GENERATE_SUPERSAMPLE_QC`** (`modules/generate_supersample_qc.nf`, wraps `bin/generate_supersample_qc.py`) — always runs; combines all per-subsample QC files (via `.collect()`) with the (possibly-`NO_FILE`) guide assignments into the final supersample report.

Every module has a matching `stub:` block used by `-stub-run` — when editing a module's real output filenames/paths, update its `stub:` block to match or stub runs will silently diverge from real behavior.

Python scripts in `bin/` are plain argparse CLIs (`--kebab-case` flags) invoked directly from process `script:` blocks (Nextflow puts `bin/` on `PATH` automatically — no need to reference the path or `chmod +x` again after cloning, permissions are already set).

## Containers

There is one unified QC container (`docker/qc/Dockerfile`, built/pushed via `docker/build_and_push.sh` to ECR) used by `GENERATE_REPORT_DATA`, `CONCATENATE`, and `GUIDE_ASSIGNMENT` — it bundles scanpy/anndata/CRISPAT plus Node/Puppeteer and a cloned `sankeymatic_local` for report visualizations. DRAGEN runs in a separate, Illumina-provided container passed via `--dragen_container`. Both container images are passed in as required params (`--qc_container`, `--dragen_container`), not hardcoded.

## ICA export

`export_pipeline_to_ica.py` imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow pipeline (reads an API key from `~/.icav2/api_key.txt`, prompts interactively for which ICA project). `nextflow_schema.json` drives the parameter form ICA renders, so keep it in sync with any new/changed/removed `params.*` in `main.nf`.
