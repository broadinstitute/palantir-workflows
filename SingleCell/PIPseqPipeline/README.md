# PIPseq QC Nextflow Pipeline

This is a pipeline for processing single-cell QC metrics from PIPseq data using Nextflow. It runs Illumina DRAGEN's scRNA pipeline per subsample, generates per-subsample and supersample-level QC reports, and (optionally) performs CRISPR guide assignment using two independent methods: CRISPAT and a purity-based method.

The pipeline is designed to run on Illumina Connected Analytics (ICA) — DRAGEN is scheduled onto ICA's FPGA-preset pods (see `pod annotation:` lines in `modules/dragen_scrna.nf`) — but it is plain Nextflow DSL2 and can run anywhere a compatible executor and the required container images are available.

There are two entrypoints, sharing the same underlying engine (`workflows/pipseq_core.nf`):
- **`main.nf`** — production entrypoint. Describes potentially many subsamples via a `--fastq_list` CSV.
- **`main_simple.nf`** — for one-off runs with a single subsample, where hand-writing a `--fastq_list` CSV is unnecessary friction. Takes flat expression/feature/hashing FASTQ lists instead; see [Simple single-subsample entrypoint](#simple-single-subsample-entrypoint).

## Overview

### Terminology

- **Subsample**: An individual sequencing unit (e.g. a well or technical replicate), identified by `RGSM` in the fastq list.
- **Supersample**: A logical grouping of subsamples that should be analyzed together (e.g. all wells from one biological sample), identified by `--supersample_id` / `--supersample_basename`.

### Pipeline Workflow

1. **Run DRAGEN scRNA** (`DRAGEN_SCRNA`): runs once per subsample, producing per-subsample metrics, barcode summary, and filtered matrix/barcodes/features files.
2. **Generate per-subsample QC** (`GENERATE_SUBSAMPLE_QC`): always runs, one invocation per subsample, regardless of whether guide assignment is enabled.
3. **Concatenate subsamples** (`CONCATENATE`): always runs; merges all subsamples' matrices into one supersample-level AnnData (`.h5ad`) and extracts a CRISPR-features-only AnnData (`.crispr.h5ad`). Handles the single-subsample case automatically.
4. **CRISPR guide assignment** (optional; runs only if `--run_guide_assignment` is `true`, the default): two independent methods run on the same concatenated CRISPR features and publish separately —
   - **`CRISPAT_GUIDE_ASSIGNMENT`**: CRISPAT's Poisson-Gaussian mixture model.
   - **`PURITY_BASED_GUIDE_ASSIGNMENT`**: a simpler purity/count-threshold heuristic (see [Purity-Based Guide Assignment](#purity-based-guide-assignment) below).
5. **Generate supersample QC** (`GENERATE_SUPERSAMPLE_QC`): always runs; combines all per-subsample QC files with CRISPAT's guide assignment results (if available) into the final supersample-level report. The purity-based assignments are not folded into this report.

## Pipeline Structure

```
SingleCell/PIPseqPipeline/
├── main.nf                          # Production entrypoint (--fastq_list)
├── main_simple.nf                   # Simple single-subsample entrypoint (flat FASTQ params)
├── workflows/
│   └── pipseq_core.nf                # Shared engine called by both entrypoints
├── nextflow.config                  # Pipeline configuration (default params, resources, reports)
├── nextflow_schema.json             # Parameter schema for main.nf (drives the ICA-rendered input form)
├── nextflow_schema_simple.json      # Parameter schema for main_simple.nf
├── modules/
│   ├── dragen_scrna.nf               # Run DRAGEN scRNA for one subsample
│   ├── generate_subsample_qc.nf      # Generate per-subsample QC metrics
│   ├── concatenate.nf                # Concatenate subsamples into a supersample AnnData
│   ├── crispat_guide_assignment.nf   # CRISPAT guide assignment
│   ├── purity_based_guide_assignment.nf  # Purity-based guide assignment
│   └── generate_supersample_qc.nf    # Generate supersample-level QC report
├── bin/
│   ├── generate_subsample_qc.py      # Per-subsample QC metrics script
│   ├── concatenate_samples.py        # Concatenation script
│   ├── run_crispat_guide_assignment.py   # CRISPAT guide assignment script
│   ├── purity_based_guide_assignment.py  # Purity-based guide assignment script
│   └── generate_supersample_qc.py    # Supersample QC report script
├── docker/                          # Dockerfile/build scripts for the qc_container image
├── stub_test/                       # Example inputs + `-stub-run` test setup
└── README.md                        # This file
```

## Quick Start

### Prerequisites

- Nextflow (>= 25.10.0 — required by the pinned `nf-schema` validation plugin; requires network access to the Nextflow plugin registry on first run)
- A DRAGEN container image (`--dragen_container`) — provided by Illumina, not built from this repo
- A QC container image (`--qc_container`) — built from `docker/qc/Dockerfile`, see `docker/SETUP.md`
- An executor/environment that can run the `container` directive (e.g. Nextflow's k8s executor on ICA, or Docker/Singularity enabled locally via your own config)

There is no bundled Nextflow profile for local/container-less execution — every process declares a `container`, so running for real requires a container-capable executor.

### Stub run (no containers required)

Every process has a `stub:` block that just touches placeholder output files, so you can validate the pipeline's wiring without DRAGEN, CRISPAT, or any container:

```bash
cd stub_test
nextflow run ../main.nf -stub-run -params-file stub_inputs/pipeline_input.json

# Simple entrypoint
nextflow run ../main_simple.nf -stub-run -params-file stub_inputs/pipeline_input_simple.json
```

See `stub_test/stub_inputs/` for additional variants (no feature library, cell hashing).

### Running the Pipeline

```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --fastq_list fastq_list.csv \
  --supersample_id "Supersample_A" \
  --supersample_basename "supersample_a" \
  --min_valid_guides 1 \
  --max_valid_guides 2 \
  --ref_tar /path/to/reference.tar \
  --annotation_file /path/to/annotation.gtf \
  --dragen_container <dragen image> \
  --qc_container <qc image> \
  --outdir results
```

All resource params (CPU/memory/scratch-space allocations per process, plus `dragen_machine_type`) have defaults in `nextflow.config` and don't need to be passed unless you want to override them — see the **Optional** list under [Command-Line Options](#command-line-options).

### fastq_list format

CSV file with columns `RGID, RGSM, RGTY, Read1File, Read2File`:

```csv
RGID,RGSM,RGTY,Read1File,Read2File
lib1_expr,Subsample_001,expression,/path/to/lib1_R1.fastq.gz,/path/to/lib1_R2.fastq.gz
lib1_feat,Subsample_001,feature,/path/to/lib1_feat_R1.fastq.gz,/path/to/lib1_feat_R2.fastq.gz
lib2_expr,Subsample_002,expression,/path/to/lib2_R1.fastq.gz,/path/to/lib2_R2.fastq.gz
```

- `RGSM` values are subsample IDs — all rows with the same `RGSM` belong to the same subsample and are passed to DRAGEN together.
- `RGTY` indicates readgroup type: `expression`, `feature` (CRISPR/feature-barcode library), or `hashing` (cell-hashing library). A subsample can mix multiple `RGTY` values across rows.

### Simple single-subsample entrypoint

For one-off runs with a single subsample, `main_simple.nf` skips the `--fastq_list` CSV entirely — you pass expression/feature/hashing FASTQs directly, and `subsample_id` is set to `--supersample_id`:

```bash
nextflow run main_simple.nf \
  --num_input_cells 10000 \
  --expression_r1_fastqs r1_lane1.fastq.gz,r1_lane2.fastq.gz \
  --expression_r2_fastqs r2_lane1.fastq.gz,r2_lane2.fastq.gz \
  --supersample_id "Sample_A" \
  --supersample_basename "sample_a" \
  --min_valid_guides 1 \
  --max_valid_guides 2 \
  --ref_tar /path/to/reference.tar \
  --annotation_file /path/to/annotation.gtf \
  --dragen_container <dragen image> \
  --qc_container <qc image> \
  --outdir results
```

- `--expression_r1_fastqs` / `--expression_r2_fastqs` are **required** lists of files, matched by position (supports multiple lanes — just list multiple files).
- `--feature_r1_fastqs` / `--feature_r2_fastqs` are optional; if given, `--scrna_feature_barcode_reference` is required.
- `--hashing_r1_fastqs` / `--hashing_r2_fastqs` are optional; if given, `--scrna_cell_hashing_reference` is required.
- Every other param (`num_input_cells`, `min_valid_guides`/`max_valid_guides`, `ref_tar`, `annotation_file`, `dragen_container`, `qc_container`, `run_guide_assignment`, resource params, etc.) is identical to `main.nf` — see [Command-Line Options](#command-line-options) below.
- Internally, `main_simple.nf` synthesizes a DRAGEN-compatible fastq-list CSV from the given FASTQ lists and hands it to the same shared engine (`workflows/pipseq_core.nf`) `main.nf` uses — everything downstream of subsample discovery (DRAGEN, concatenation, guide assignment, QC) behaves identically either way.
- Validated against `nextflow_schema_simple.json` (a separate schema from `main.nf`'s `nextflow_schema.json`, since the input params differ).

### Parameter validation

Required/typed params (presence, type, allowed range/pattern) are validated against `nextflow_schema.json` via the [`nf-schema`](https://nextflow-io.github.io/nf-schema/) plugin as soon as the pipeline starts — a missing, mistyped, or out-of-range param fails immediately with a clear message rather than partway through the run.

A few checks that can't be expressed in JSON Schema are enforced separately, right after schema validation:
- `RGTY` values in `--fastq_list` must be exactly `expression`, `feature`, or `hashing` (case-sensitive) — a typo fails immediately instead of silently producing an empty feature/hashing group.
- If `--fastq_list` has any `feature` rows, `--scrna_feature_barcode_reference` must be set (and likewise `--scrna_cell_hashing_reference` for `hashing` rows).
- `--min_valid_guides` must be `<= --max_valid_guides`.

### Command-Line Options

For `main.nf` (`main_simple.nf` shares everything here except `--fastq_list`, which it replaces with the FASTQ-list params described in [Simple single-subsample entrypoint](#simple-single-subsample-entrypoint)):

**Required:**
- `--num_input_cells`: Number of input cells (integer)
- `--fastq_list`: CSV file described above
- `--supersample_id`: Supersample identifier
- `--supersample_basename`: Supersample basename for output organization
- `--min_valid_guides` / `--max_valid_guides`: Guide-count thresholds used for guide assignment QC (integers; `0` is a valid value for `--min_valid_guides`)
- `--ref_tar`: DRAGEN reference genome tar file
- `--annotation_file`: Gene annotation file (GTF/GFF) for DRAGEN
- `--dragen_container`: Container image for DRAGEN execution
- `--qc_container`: Container image for QC processing

**Optional:**
- `--run_guide_assignment`: Whether to run CRISPR guide assignment (default: `true`)
- `--use_direct_capture_mode`: Whether to use DRAGEN direct-capture mode for feature barcodes (default: `true`)
- `--scrna_feature_barcode_reference`: Feature barcode reference CSV for DRAGEN (only needed if the fastq_list has `feature` rows)
- `--scrna_barcode_sequence_list`: Barcode sequence list CSV for DRAGEN
- `--scrna_cell_hashing_reference`: Cell hashing reference CSV for DRAGEN (only needed if the fastq_list has `hashing` rows)
- `--additional_dragen_args`: Extra raw arguments appended to the DRAGEN command line
- `--guide_assignment_num_processes`: Number of processes for CRISPAT guide assignment (default: all available cores)
- `--concatenate_cpus` / `--concatenate_memory_gb`: Resources for the `CONCATENATE` process (defaults: `16` / `64`)
- `--dragen_scratch_tb`: Scratch disk space (TiB) allocated to the DRAGEN pod (default: `2`)
- `--dragen_machine_type`: ICA pod preset size for the DRAGEN job (default: `fpga2-medium`)
- `--cpu_generate_subsample_qc` / `--memory_gb_generate_subsample_qc`: Resources for the `GENERATE_SUBSAMPLE_QC` process (defaults: `8` / `32`)
- `--cpu_generate_supersample_qc` / `--memory_gb_generate_supersample_qc`: Resources for the `GENERATE_SUPERSAMPLE_QC` process (defaults: `8` / `32`)
- `--cpu_crispat_guide_assignment` / `--memory_gb_crispat_guide_assignment`: Resources for the `CRISPAT_GUIDE_ASSIGNMENT` process (defaults: `8` / `32`)
- `--cpu_purity_based_guide_assignment` / `--memory_gb_purity_based_guide_assignment`: Resources for the `PURITY_BASED_GUIDE_ASSIGNMENT` process (defaults: `8` / `32`)
- `--outdir`: Output directory (default: `out`)
- `--help`: Show help message

## Output

Results are organized under `${params.outdir}/${params.supersample_basename}/`:

- **`<subsample_id>/dragen_output/`**: raw DRAGEN outputs for that subsample (metrics CSV, barcode summary, filtered matrix/barcodes/features)
- **`<subsample_id>/logs/`**: DRAGEN logs for that subsample
- **`<subsample_id>/qc/`**: per-subsample QC files
  - `<subsample_id>.qc_metrics.tsv`
  - `<subsample_id>.qc_barcode_metrics.tsv`
- **`adata/`**: concatenated AnnData files (always produced, independent of guide assignment)
  - `<supersample_basename>.h5ad` — full concatenated dataset
  - `<supersample_basename>.crispr.h5ad` — CRISPR-features-only subset
- **`crispat_ga/`**: CRISPAT guide assignment outputs (only if `--run_guide_assignment true`)
  - `poisson_gauss/assignments.csv`
- **`purity_ga/`**: purity-based guide assignment output (only if `--run_guide_assignment true`)
  - `<supersample_id>.purity_based_guide_assignments.csv`
- **`supersample_qc/`**: final supersample-level report
  - `<supersample_basename>.supersample_qc_metrics.tsv`
  - `<supersample_basename>.guide_assignment_distribution.png` (only if guide assignment ran)
- **`pipeline_info/`**: Nextflow reports (`timeline.html`, `report.html`, `trace.txt`, `dag.svg`)

A `README.txt` describing this layout is written directly into `${params.outdir}/${params.supersample_basename}/` when the run finishes successfully.

## CRISPR Guide Assignment

Guide assignment is **enabled by default** (`--run_guide_assignment true`) and adds two independent steps on top of the concatenation that always happens, both consuming the same `<supersample_basename>.crispr.h5ad`:

1. **Concatenate**: merge all subsamples and extract CRISPR Direct Capture features into `<supersample_basename>.crispr.h5ad` (`bin/concatenate_samples.py`) — this always runs.
2. **CRISPAT guide assignment**: run CRISPAT's Poisson-Gaussian mixture model on the CRISPR AnnData (`bin/run_crispat_guide_assignment.py`) — only if enabled. Feeds into the final supersample QC report.
3. **Purity-based guide assignment**: an independent, simpler heuristic (`bin/purity_based_guide_assignment.py`) — only if enabled. Published on its own; not folded into the supersample QC report.
4. **Supersample report**: fold CRISPAT's guide assignment results into the final QC report (`bin/generate_supersample_qc.py`) — always runs, with or without guide assignment data.

Disable both guide assignment methods with:
```bash
nextflow run main.nf ... --run_guide_assignment false
```

### Purity-Based Guide Assignment

For each cell, `bin/purity_based_guide_assignment.py` looks at the CRISPR guide UMI counts and computes:
- `total_count`: sum of all guide counts for that cell
- `count_1st` / `count_2nd`: the highest and second-highest guide counts for that cell
- `purity`: `count_1st / (count_1st + count_2nd)`

A cell is assigned to its top guide (`gRNA`) only if `total_count > 10` **and** `purity > 0.75`; otherwise `gRNA` is left empty. Output columns: `cell, gRNA, purity, total_count, count_1st, count_2nd`.

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run main.nf ... -resume
```

## Troubleshooting

**Error**: `ModuleNotFoundError` / missing Python package inside a process
- **Solution**: The package needs to be added to `docker/qc/Dockerfile` (or `docker/qc/requirements.txt`) and the `qc_container` image rebuilt/pushed — pipeline processes run inside the container, not your local environment.

**Error**: CRISPR features extraction fails / no CRISPR features found
- **Solution**: Verify your `scrna_feature_barcode_reference` and fastq_list `feature` rows are correct, and that DRAGEN's features output actually contains a "CRISPR Direct Capture" feature type.

**Error**: Pipeline exits immediately with a parameter validation error
- **Solution**: Required/typed params are validated against `nextflow_schema.json` via the `nf-schema` plugin (see [Parameter validation](#parameter-validation)); a few additional business-logic checks (RGTY values, guide-count range) run right after. Check the exact list of required flags above.

**Out of memory / timeout**: Adjust resource allocations in `nextflow.config` for specific processes, e.g.:
  ```groovy
  process {
      withName: CRISPAT_GUIDE_ASSIGNMENT {
          memory = 64.GB
          time = 24.h
      }
  }
  ```

## Additional Resources

- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
- [Nextflow Container Documentation](https://www.nextflow.io/docs/latest/container.html)
- [CRISPAT Documentation](https://github.com/pinellolab/CRISPAT)
