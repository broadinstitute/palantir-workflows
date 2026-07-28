# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow DSL2 pipeline for processing PIPseq single-cell RNA-seq data: it runs DRAGEN scRNA per subsample, then produces per-subsample and supersample-level QC reports, with optional CRISPR guide assignment via two independent methods (CRISPAT and a purity-based heuristic). It's part of the larger `palantir-workflows` monorepo (see `../../.github/copilot-instructions.md` for repo-wide conventions) but is developed and run from this directory. The deployment target is Illumina Connected Analytics (ICA) — DRAGEN runs on ICA's FPGA-preset pods (see `pod annotation:` lines in `modules/dragen_scrna.nf`), and `ica_tools/export_pipeline_to_ica.py` (plus the `nextflow_schema*.json` files, which drive the ICA-rendered input forms) exist to publish this pipeline into ICA.

There are **two entrypoints**, sharing one engine:
- **`main.nf`** — production entrypoint, describes potentially many subsamples via a `--fastq_list` CSV. Validated against `nextflow_schema.json`.
- **`main_simple.nf`** — one-off single-subsample entrypoint, takes flat expression/feature/hashing FASTQ list params instead of a CSV, and sets `subsample_id = supersample_id`. Validated against `nextflow_schema_simple.json`.
- Both call into `workflows/pipseq_core.nf`'s named `PIPSEQ_CORE` workflow, which contains everything from DRAGEN onward — see [Architecture](#architecture).

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
nextflow run ../main_simple.nf -stub-run -params-file stub_inputs/pipeline_input_simple.json

# Real run (main.nf)
nextflow run main.nf --num_input_cells <int> --fastq_list <fastq_list.csv> \
  --supersample_id <id> --supersample_basename <name> \
  --min_valid_guides <int> --max_valid_guides <int> \
  --ref_tar <ref.tar> --annotation_file <annotation.gtf> \
  --dragen_container <image> --qc_container <image> \
  --concatenate_cpus <n> --concatenate_memory_gb <n> --dragen_scratch_tb <n>
```

`stub_test/stub_inputs/` has several ready-made `pipeline_input*.json` variants (with/without feature barcodes, with cell hashing, and the simple entrypoint) for exercising different branches of the parsing logic.

## Parameter validation

Each entrypoint has its own schema (`main.nf` → `nextflow_schema.json`, `main_simple.nf` → `nextflow_schema_simple.json`, passed via `validateParameters(parameters_schema: '...')`) — these are the single source of truth for "what's required" for that entrypoint. Both `include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'` and call `validateParameters()` right after their own `--help` check, which validates presence/type/pattern/min-max for every param against the schema and fails fast with a clear message. **Don't add hand-rolled `if (!params.x) { exit 1 }` checks for anything expressible in JSON Schema** — add/edit the corresponding property (and its `required` list) in the relevant schema file instead, so ICA's rendered form and the pipeline's own validation never drift apart again (they did, historically). The two schemas share almost all of their non-`input_output_options` sections (`dragen_options`, `guide_assignment_options`, `container_options`, `resource_options`, `generic_options`) verbatim — keep both in sync when changing a shared param.

A few checks that JSON Schema can't express are hand-written per entrypoint, right after `validateParameters()`:
- `main.nf`: `RGTY` values in the fastq_list must be one of `expression`/`feature`/`hashing` (checked once the CSV is fully parsed, inside the `subsample_info` `.flatMap` block); `--scrna_feature_barcode_reference` must be set if any row has `RGTY=feature` (same for `--scrna_cell_hashing_reference` / `RGTY=hashing`).
- `main_simple.nf`: each FASTQ category's R1/R2 list lengths must match; the same feature/hashing-reference cross-checks, phrased as "were feature/hashing FASTQs given" instead of "does the CSV have feature/hashing rows".
- Both, inside `workflows/pipseq_core.nf` (genuinely entrypoint-agnostic): `--min_valid_guides <= --max_valid_guides`.

**Nextflow `params` gotcha:** session-level params (set via CLI/`-params-file`/schema defaults at launch) are visible identically from every included module. But a plain `params.x = someValue` assignment made at runtime, inside an entrypoint's `workflow {}` body, does **not** propagate across an `include` boundary — `workflows/pipseq_core.nf` would still see the old/default value. This is why the fastq-list *path* is threaded into `PIPSEQ_CORE` as an explicit `take:` parameter (`fastq_list_path`) rather than via `params.fastq_list` — `main_simple.nf` synthesizes a CSV at runtime and has no other way to hand that path to the shared workflow. If you need to pass an entrypoint-computed value into `PIPSEQ_CORE`, add it to `take:`; don't rely on a runtime `params.x = ...` write being visible there.

## Architecture

Both entrypoints are thin: parse their own input shape into a `subsample_info` list (one map per subsample: `[rgsm, feature_rgids, hashing_rgids, fastq_files]`), plus (for `main_simple.nf`) synthesize a DRAGEN-compatible fastq-list CSV, then call `PIPSEQ_CORE(subsample_info, fastq_list_path)` from `workflows/pipseq_core.nf`. Everything from DRAGEN onward lives in that one shared named workflow and doesn't know or care which entrypoint built `subsample_info`. Data flows through modules in `modules/`:

0. **Entrypoint-specific subsample discovery:**
   - `main.nf`: **Parse `--fastq_list`** (CSV with `RGID, RGSM, RGTY, Read1File, Read2File`). Rows are grouped by `RGSM` (subsample). `RGTY` is `'expression'`, `'feature'`, or `'hashing'`; feature/hashing RGIDs are collected per-subsample into comma-joined lists (`feature_barcode_groups`, `hto_barcode_groups`). `--fastq_files` (required, unused by any pipeline logic) exists purely so ICA localizes the FASTQ files referenced by path *inside* the `--fastq_list` CSV onto the compute node before the run starts — ICA only localizes files declared as their own top-level input, and has no way to know the CSV references other files. It's required rather than optional on the assumption that this pipeline currently only ever runs on ICA. `main_simple.nf` doesn't need an equivalent, since its FASTQ params (`expression_r1_fastqs` etc.) are already individually-declared top-level inputs.
   - `main_simple.nf`: builds one synthetic row per (R1, R2) pair across the `expression`/`feature`/`hashing` FASTQ list params, all sharing `RGSM = supersample_id` (always exactly one subsample), writes them out as a DRAGEN-compatible CSV (`RGID,RGSM,RGLB,Lane,Read1File,Read2File,RGTY` — Illumina's `--fastq-list` spec requires `RGID,RGSM,RGLB,Lane,Read1File,Read2File`; any other 4-uppercase-char `RG**`-named column, like `RGTY`, is a supported passthrough custom tag) into a file under `workflow.workDir`, and builds the matching single-element `subsample_info` directly from the same file lists (no need to re-parse the CSV it just wrote).
1. **`PIPSEQ_CORE`** (`workflows/pipseq_core.nf`) — the shared engine both entrypoints call. Feature/hashing RGIDs from `subsample_info` get bundled into a per-subsample **meta map** (see below) passed to `DRAGEN_SCRNA`.
2. **`DRAGEN_SCRNA`** (`modules/dragen_scrna.nf`) — one invocation per subsample. Declares two kinds of outputs: a full `path 'dragen_output/*', emit: output` glob that publishes **every** file DRAGEN writes (this must stay — it's what guarantees nothing DRAGEN produces is silently dropped just because this pipeline doesn't read it), plus named per-file-type emits (`metrics`, `barcode_summary`, `matrix`, `barcodes`, `features`), each a `(subsample_id, file)` tuple. `PIPSEQ_CORE` rejoins those named channels with `.join()` keyed on `subsample_id` — there's no filename parsing anymore; if you rename one of the 5 filenames DRAGEN produces, update the corresponding `path("dragen_output/${meta.subsample_id}...")` output declaration (and its `stub:` block) rather than any string-matching logic.
3. **`GENERATE_SUBSAMPLE_QC`** (`modules/generate_subsample_qc.nf`, wraps `bin/generate_subsample_qc.py`) — runs per subsample, always executes regardless of guide assignment.
4. **`CONCATENATE`** (`modules/concatenate.nf`, wraps `bin/concatenate_samples.py`) — always runs (independent of `run_guide_assignment`); merges all subsamples' matrices into one supersample-level `.h5ad` plus a CRISPR-features-only `.crispr.h5ad`. Handles the single-subsample case automatically (no special-casing needed by the caller).
5. **Guide assignment** — only runs if `params.run_guide_assignment` (default `true`); two independent methods both consume `CONCATENATE`'s CRISPR AnnData output and publish separately, neither feeding the other:
   - **`CRISPAT_GUIDE_ASSIGNMENT`** (`modules/crispat_guide_assignment.nf`, wraps `bin/run_crispat_guide_assignment.py`) — CRISPAT's Poisson-Gaussian mixture model. Its output (`crispat_guide_assignments_ch` in `workflows/pipseq_core.nf`) is what feeds `GENERATE_SUPERSAMPLE_QC` below.
   - **`PURITY_BASED_GUIDE_ASSIGNMENT`** (`modules/purity_based_guide_assignment.nf`, wraps `bin/purity_based_guide_assignment.py`) — a vectorized numpy/scanpy heuristic, no CRISPAT dependency. For each cell: `total_count` (sum of all guide counts), `count_1st`/`count_2nd` (top two guide counts), `purity = count_1st / (count_1st + count_2nd)`. Assigns the top guide as `gRNA` only if `total_count > 10` and `purity > 0.75`; otherwise `gRNA` is empty. Not consumed by any other process — published on its own under `purity_ga/`.

   When `run_guide_assignment` is disabled, downstream code substitutes a `Channel.of(file('NO_FILE'))` sentinel for `crispat_guide_assignments_ch` — modules check `guide_assignments.name != 'NO_FILE'` rather than checking optionality via Nextflow's `optional` input mechanism.
6. **`GENERATE_SUPERSAMPLE_QC`** (`modules/generate_supersample_qc.nf`, wraps `bin/generate_supersample_qc.py`) — always runs; combines all per-subsample QC files (via `.collect()`) with the (possibly-`NO_FILE`) CRISPAT guide assignments into the final supersample report.

On successful completion, `workflow.onComplete` writes a plain-text `README.txt` into the output directory describing this layout.

Every module has a matching `stub:` block used by `-stub-run` — when editing a module's real output filenames/paths, update its `stub:` block to match or stub runs will silently diverge from real behavior.

Python scripts in `bin/` are plain argparse CLIs (`--kebab-case` flags) invoked directly from process `script:` blocks (Nextflow puts `bin/` on `PATH` automatically — no need to reference the path or `chmod +x` again after cloning, permissions are already set).

### Meta-map convention

`DRAGEN_SCRNA`, `GENERATE_SUBSAMPLE_QC`, and `GENERATE_SUPERSAMPLE_QC` take a `tuple val(meta), path(...), ...` input, following the common nf-core convention: scalar/string values that don't need file staging are bundled into a `meta` Groovy map (e.g. `meta.subsample_id`, `meta.num_input_cells`), while anything that needs to be staged into the task's work directory stays as an explicit named `path()` element. **Don't hide `Path`/`file()` objects inside a meta map that gets passed into a process call** — Nextflow only stages files declared via `path()`; a `Path` sitting inside a plain `val(map)` won't be staged and will break under non-local executors (S3/GS/remote work dirs, containers). It's fine to bundle `Path`s into a Groovy map purely as an in-memory convenience *between* channel operators in an entrypoint or in `workflows/pipseq_core.nf`, as long as the map gets unpacked back into individual `path()`/`val()` tuple elements before being handed to the next process call. `CONCATENATE`, `CRISPAT_GUIDE_ASSIGNMENT`, and `PURITY_BASED_GUIDE_ASSIGNMENT` don't use this pattern since they're aggregate steps over all subsamples (or take a single simple `path` input), not naturally "one meta per call".

`modules/dragen_scrna.nf` also defines a module-level `buildDragenOptionalArgs()` Groovy function, shared between `script:` and `stub:`, so the conditional-argument logic for optional DRAGEN inputs (feature/hashing/barcode-list) can't drift between the two.

## Containers

There is one unified QC container (`docker/qc/Dockerfile`, built/pushed via `docker/build_and_push.sh` to ECR) used by `GENERATE_SUBSAMPLE_QC`, `CONCATENATE`, `CRISPAT_GUIDE_ASSIGNMENT`, `PURITY_BASED_GUIDE_ASSIGNMENT`, and `GENERATE_SUPERSAMPLE_QC` — it bundles scanpy/anndata/CRISPAT plus Node/Puppeteer and a cloned `sankeymatic_local` for a currently-disabled Sankey-plot feature (see the commented-out `generate_sankey_plot()` call and `N2` computation in `bin/generate_supersample_qc.py` — left in place intentionally in case that feature gets finished later). DRAGEN runs in a separate, Illumina-provided container passed via `--dragen_container`. Both container images are passed in as required params (`--qc_container`, `--dragen_container`), not hardcoded.

## ICA export

All ICA-facing tooling lives in `ica_tools/`: `export_pipeline_to_ica.py`, `start_analysis.py`, `ica_common.py` (shared helpers), and `ica_tools/inputforms/<main|main_simple>/inputForm.json` (hand-maintained ICA launch-form definitions, one per entrypoint).

`ica_tools/export_pipeline_to_ica.py` imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow pipeline (reads an API key from `~/.icav2/api_key.txt`, prompts interactively for which ICA project **and which entrypoint** — `main.nf` or `main_simple.nf` — since they're registered as separate ICA pipelines sharing the same `nextflow.config`). Each `nextflow_schema*.json` drives the parameter form ICA renders for its entrypoint, so keep it in sync with any new/changed/removed `params.*` in the corresponding entrypoint file (or in `workflows/pipseq_core.nf`, for shared params) — this is now enforced at pipeline-runtime too via `validateParameters()`, not just relevant to ICA's form.

After the `:importGitPipeline` call, the script polls `GET /projects/{projectId}/pipelines/{pipelineId}` (via `pipeline.statusAsString`) until the git import finishes and the pipeline reaches `Draft` status, then `PUT`s the corresponding `ica_tools/inputforms/<main|main_simple>/inputForm.json` to `/projects/{projectId}/pipelines/{pipelineId}/inputForm/inputFormFile` — this is what previously had to be pasted in by hand via the ICA web UI's pipeline editor. **Keep `ica_tools/inputforms/<entrypoint>/inputForm.json` in sync whenever a param is added/changed/removed** in the corresponding `nextflow_schema*.json` — it's a separately hand-maintained file (ICA doesn't derive its launch form from the schema), so nothing enforces this at runtime the way `validateParameters()` does for the schema itself.

`ica_common.py` holds what both ICA scripts share: the API URL, the `~/.icav2/api_key.txt` reader, the (`BCL Shared Development` / `MDL Single Cell Dev`) project-ID map, and the `prompt_choice()` numbered-menu helper.

`ica_tools/start_analysis.py` (no CLI args other than `--dry-run`) interactively starts an actual analysis run for a pipeline already imported into ICA: pick the ICA project, pick which already-imported PIPseq pipeline to run (listed via `GET /projects/{projectId}/pipelines`, filtered to pipelines whose `code` starts with `PIPseq_BCL`, newest first, flagging whichever one's `gitPipelineImportDto.commitId` matches the current git HEAD), then pick `test/test_inputs_main.json` or `test/test_inputs_main_simple.json` — each is *just* the flat `inputs` map for that entrypoint (same shape as `stub_test/stub_inputs/pipeline_input*.json`), no project/pipeline metadata. `--dry-run` resolves inputs and prints the request payload without submitting. File/folder inputs are given as project-relative paths (e.g. `/some_folder/some_file.fastq.gz`) rather than ICA data IDs — the script resolves each path to its `fil.<hash>`/`fol.<hash>` ID via `GET /projects/{projectId}/data?filePath=...&filePathMatchMode=FULL_CASE_INSENSITIVE`, using the *chosen pipeline's own live* input form (fetched from `GET .../pipelines/{pipelineId}/inputForm/inputFormFile` — the same endpoint `export_pipeline_to_ica.py` uploads to) to know which fields are `"type": "data"`, rather than the local `ica_tools/inputforms/*.json` copy, which could drift from what's actually deployed (a value that already looks like a `fil./fol.` ID is passed through unresolved). Submits via `POST /projects/{projectId}/analysis:nextflowJson`. Note some ICA endpoints only accept specific `application/vnd.illumina.v*+json`/`application/octet-stream` Accept values (`analysis:nextflowJson` wants v4 JSON, the data-listing endpoint wants v3 JSON, the input-form download wants octet-stream) — mismatches fail with a 400 "Invalid Accept Header" (see the `inputForm/inputFormFile` upload bug this bit us with earlier).
