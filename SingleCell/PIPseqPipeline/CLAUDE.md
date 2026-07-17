# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow DSL2 pipeline for processing PIPseq single-cell RNA-seq data: it runs DRAGEN scRNA per subsample, then produces per-subsample and supersample-level QC reports, with optional CRISPR guide assignment via two independent methods (CRISPAT and a purity-based heuristic). It's part of the larger `palantir-workflows` monorepo (see `../../.github/copilot-instructions.md` for repo-wide conventions) but is developed and run from this directory. The deployment target is Illumina Connected Analytics (ICA) — DRAGEN runs on ICA's FPGA-preset pods (see `pod annotation:` lines in `modules/dragen_scrna.nf`), and `export_pipeline_to_ica.py` (plus `nextflow_schema.json`, which drives the ICA-rendered input form) exist to publish this pipeline into ICA.

`README.md` and `docker/README.md` describe the actual current CLI/params/behavior — trust them (and `main.nf`/`nextflow_schema.json`/`modules/*.nf`) over any older notes you find elsewhere in the repo (e.g. `doc/`).

## Terminology

- **Subsample**: an individual sequencing unit (e.g. one well/technical replicate), identified by `RGSM` in the fastq list.
- **Supersample**: a logical grouping of subsamples analyzed together, identified by `--supersample_id` / `--supersample_basename`.

## Running / testing

There's no unit test suite; correctness is verified by running the pipeline (or Nextflow's `-stub-run` mode) end to end. Nextflow **>=25.10.0** is required (pinned by the `nf-schema` plugin in `nextflow.config`); first run needs network access to the Nextflow plugin registry to download it.

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

`stub_test/stub_inputs/` has several ready-made `pipeline_input*.json` variants (with/without feature barcodes, with cell hashing) for exercising different branches of the fastq-list parsing logic.

## Parameter validation

`nextflow_schema.json` is the single source of truth for "what's required" — `main.nf` calls `validateParameters()` (from the `nf-schema` plugin, included via `include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'`) right after the `--help` check, which validates presence/type/pattern/min-max for every param against the schema and fails fast with a clear message. **Don't add hand-rolled `if (!params.x) { exit 1 }` checks for anything expressible in JSON Schema** — add/edit the corresponding property (and its `required` list) in `nextflow_schema.json` instead, so ICA's rendered form and the pipeline's own validation never drift apart again (they did, historically).

A few checks that JSON Schema can't express are hand-written in `main.nf` right after `validateParameters()`:
- `RGTY` values in the fastq_list must be one of `expression`/`feature`/`hashing` (checked once the CSV is fully parsed, inside the `subsample_info` `.flatMap` block).
- `--scrna_feature_barcode_reference` must be set if any row has `RGTY=feature` (same for `--scrna_cell_hashing_reference` / `RGTY=hashing`).
- `--min_valid_guides <= --max_valid_guides`.
- `--scrna_feature_barcode_reference` can't contain more than `MAX_CRISPR_GUIDES` (300) guides — checked before any DRAGEN job runs. `bin/concatenate_samples.py` re-checks the same limit after concatenation as a safety net (in case DRAGEN's actual CRISPR feature count ever disagrees with the reference file's guide count).

## Architecture

`main.nf` is a single linear `workflow {}` block (no subworkflows). Data flows through modules in `modules/`:

1. **Parse `--fastq_list`** (CSV with `RGID, RGSM, RGTY, Read1File, Read2File`). Rows are grouped by `RGSM` (subsample). `RGTY` is `'expression'`, `'feature'`, or `'hashing'`; feature/hashing RGIDs are collected per-subsample into comma-joined lists (`feature_barcode_groups`, `hto_barcode_groups`) that get bundled into a per-subsample **meta map** (see below) passed to `DRAGEN_SCRNA`.
2. **`DRAGEN_SCRNA`** (`modules/dragen_scrna.nf`) — one invocation per subsample. Declares two kinds of outputs: a full `path 'dragen_output/*', emit: output` glob that publishes **every** file DRAGEN writes (this must stay — it's what guarantees nothing DRAGEN produces is silently dropped just because this pipeline doesn't read it), plus named per-file-type emits (`metrics`, `barcode_summary`, `matrix`, `barcodes`, `features`), each a `(subsample_id, file)` tuple. `main.nf` rejoins those named channels with `.join()` keyed on `subsample_id` — there's no filename parsing anymore; if you rename one of the 5 filenames DRAGEN produces, update the corresponding `path("dragen_output/${meta.subsample_id}...")` output declaration (and its `stub:` block) rather than any string-matching logic.
3. **`GENERATE_REPORT_DATA`** (`modules/generate_report_data.nf`, wraps `bin/generate_report_data.py`) — runs per subsample, always executes regardless of guide assignment.
4. **`CONCATENATE`** (`modules/concatenate.nf`, wraps `bin/concatenate_samples.py`) — always runs (independent of `run_guide_assignment`); merges all subsamples' matrices into one supersample-level `.h5ad` plus a CRISPR-features-only `.crispr.h5ad`. Handles the single-subsample case automatically (no special-casing needed by the caller).
5. **Guide assignment** — only runs if `params.run_guide_assignment` (default `true`); two independent methods both consume `CONCATENATE`'s CRISPR AnnData output and publish separately, neither feeding the other:
   - **`CRISPAT_GUIDE_ASSIGNMENT`** (`modules/crispat_guide_assignment.nf`, wraps `bin/run_crispat_guide_assignment.py`) — CRISPAT's Poisson-Gaussian mixture model. Its output (`crispat_guide_assignments_ch` in `main.nf`) is what feeds `GENERATE_SUPERSAMPLE_QC` below.
   - **`PURITY_BASED_GUIDE_ASSIGNMENT`** (`modules/purity_based_guide_assignment.nf`, wraps `bin/purity_based_guide_assignment.py`) — a vectorized numpy/scanpy heuristic, no CRISPAT dependency. For each cell: `total_count` (sum of all guide counts), `count_1st`/`count_2nd` (top two guide counts), `purity = count_1st / (count_1st + count_2nd)`. Assigns the top guide as `gRNA` only if `total_count > 10` and `purity > 0.75`; otherwise `gRNA` is empty. Not consumed by any other process — published on its own under `purity_ga/`.

   When `run_guide_assignment` is disabled, downstream code substitutes a `Channel.of(file('NO_FILE'))` sentinel for `crispat_guide_assignments_ch` — modules check `guide_assignments.name != 'NO_FILE'` rather than checking optionality via Nextflow's `optional` input mechanism.
6. **`GENERATE_SUPERSAMPLE_QC`** (`modules/generate_supersample_qc.nf`, wraps `bin/generate_supersample_qc.py`) — always runs; combines all per-subsample QC files (via `.collect()`) with the (possibly-`NO_FILE`) CRISPAT guide assignments into the final supersample report.

On successful completion, `workflow.onComplete` writes a plain-text `README.txt` into the output directory describing this layout.

Every module has a matching `stub:` block used by `-stub-run` — when editing a module's real output filenames/paths, update its `stub:` block to match or stub runs will silently diverge from real behavior.

Python scripts in `bin/` are plain argparse CLIs (`--kebab-case` flags) invoked directly from process `script:` blocks (Nextflow puts `bin/` on `PATH` automatically — no need to reference the path or `chmod +x` again after cloning, permissions are already set).

### Meta-map convention

`DRAGEN_SCRNA`, `GENERATE_REPORT_DATA`, and `GENERATE_SUPERSAMPLE_QC` take a `tuple val(meta), path(...), ...` input, following the common nf-core convention: scalar/string values that don't need file staging are bundled into a `meta` Groovy map (e.g. `meta.subsample_id`, `meta.num_input_cells`), while anything that needs to be staged into the task's work directory stays as an explicit named `path()` element. **Don't hide `Path`/`file()` objects inside a meta map that gets passed into a process call** — Nextflow only stages files declared via `path()`; a `Path` sitting inside a plain `val(map)` won't be staged and will break under non-local executors (S3/GS/remote work dirs, containers). It's fine to bundle `Path`s into a Groovy map purely as an in-memory convenience *between* channel operators in `main.nf`, as long as the map gets unpacked back into individual `path()`/`val()` tuple elements before being handed to the next process call. `CONCATENATE`, `CRISPAT_GUIDE_ASSIGNMENT`, and `PURITY_BASED_GUIDE_ASSIGNMENT` don't use this pattern since they're aggregate steps over all subsamples (or take a single simple `path` input), not naturally "one meta per call".

`modules/dragen_scrna.nf` also defines a module-level `buildDragenOptionalArgs()` Groovy function, shared between `script:` and `stub:`, so the conditional-argument logic for optional DRAGEN inputs (feature/hashing/barcode-list) can't drift between the two.

## Containers

There is one unified QC container (`docker/qc/Dockerfile`, built/pushed via `docker/build_and_push.sh` to ECR) used by `GENERATE_REPORT_DATA`, `CONCATENATE`, `CRISPAT_GUIDE_ASSIGNMENT`, `PURITY_BASED_GUIDE_ASSIGNMENT`, and `GENERATE_SUPERSAMPLE_QC` — it bundles scanpy/anndata/CRISPAT plus Node/Puppeteer and a cloned `sankeymatic_local` for a currently-disabled Sankey-plot feature (see the commented-out `generate_sankey_plot()` call and `N2` computation in `bin/generate_supersample_qc.py` — left in place intentionally in case that feature gets finished later). DRAGEN runs in a separate, Illumina-provided container passed via `--dragen_container`. Both container images are passed in as required params (`--qc_container`, `--dragen_container`), not hardcoded.

## ICA export

`export_pipeline_to_ica.py` imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow pipeline (reads an API key from `~/.icav2/api_key.txt`, prompts interactively for which ICA project). `nextflow_schema.json` drives the parameter form ICA renders, so keep it in sync with any new/changed/removed `params.*` in `main.nf` — this is now enforced at pipeline-runtime too via `validateParameters()`, not just relevant to ICA's form.
