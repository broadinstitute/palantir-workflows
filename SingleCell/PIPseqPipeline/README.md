# PIPseq QC Nextflow Pipeline

This is a pipeline for processing single-cell QC metrics from PIPseq data using Nextflow. It runs Illumina DRAGEN's scRNA pipeline per subsample, generates per-subsample and supersample-level QC reports, and (optionally) performs CRISPR guide assignment via CRISPAT.

The pipeline is designed to run on Illumina Connected Analytics (ICA) — DRAGEN is scheduled onto ICA's FPGA-preset pods (see `pod annotation:` lines in `modules/dragen_scrna.nf`) — but it is plain Nextflow DSL2 and can run anywhere a compatible executor and the required container images are available.

## Overview

### Terminology

- **Subsample**: An individual sequencing unit (e.g. a well or technical replicate), identified by `RGSM` in the fastq list.
- **Supersample**: A logical grouping of subsamples that should be analyzed together (e.g. all wells from one biological sample), identified by `--supersample_id` / `--supersample_basename`.

### Pipeline Workflow

1. **Run DRAGEN scRNA** (`DRAGEN_SCRNA`): runs once per subsample, producing per-subsample metrics, barcode summary, and filtered matrix/barcodes/features files.
2. **Generate per-subsample QC** (`GENERATE_REPORT_DATA`): always runs, one invocation per subsample, regardless of whether guide assignment is enabled.
3. **Concatenate subsamples** (`CONCATENATE`): always runs; merges all subsamples' matrices into one supersample-level AnnData (`.h5ad`) and extracts a CRISPR-features-only AnnData (`.crispr.h5ad`). Handles the single-subsample case automatically.
4. **CRISPR guide assignment** (`GUIDE_ASSIGNMENT`, optional): performs CRISPAT guide assignment (Poisson-Gaussian mixture model) on the concatenated CRISPR features. Runs only if `--run_guide_assignment` is `true` (the default).
5. **Generate supersample QC** (`GENERATE_SUPERSAMPLE_QC`): always runs; combines all per-subsample QC files with guide assignment results (if available) into the final supersample-level report.

## Pipeline Structure

```
SingleCell/PIPseqPipeline/
├── main.nf                          # Main Nextflow pipeline (entrypoint)
├── nextflow.config                  # Pipeline configuration (default params, resources, reports)
├── nextflow_schema.json             # Parameter schema (drives the ICA-rendered input form)
├── modules/
│   ├── dragen_scrna.nf               # Run DRAGEN scRNA for one subsample
│   ├── generate_report_data.nf       # Generate per-subsample QC metrics
│   ├── concatenate.nf                # Concatenate subsamples into a supersample AnnData
│   ├── guide_assignment.nf           # CRISPAT guide assignment
│   └── generate_supersample_qc.nf    # Generate supersample-level QC report
├── bin/
│   ├── generate_report_data.py       # Per-subsample QC metrics script
│   ├── concatenate_samples.py        # Concatenation script
│   ├── run_guide_assignment.py       # CRISPAT guide assignment script
│   └── generate_supersample_qc.py    # Supersample QC report script
├── docker/                          # Dockerfile/build scripts for the qc_container image
├── stub_test/                       # Example inputs + `-stub-run` test setup
└── README.md                        # This file
```

## Quick Start

### Prerequisites

- Nextflow (>= 22.10.0)
- A DRAGEN container image (`--dragen_container`) — provided by Illumina, not built from this repo
- A QC container image (`--qc_container`) — built from `docker/qc/Dockerfile`, see `docker/SETUP.md`
- An executor/environment that can run the `container` directive (e.g. Nextflow's k8s executor on ICA, or Docker/Singularity enabled locally via your own config)

There is no bundled Nextflow profile for local/container-less execution — every process declares a `container`, so running for real requires a container-capable executor.

### Stub run (no containers required)

Every process has a `stub:` block that just touches placeholder output files, so you can validate the pipeline's wiring without DRAGEN, CRISPAT, or any container:

```bash
cd stub_test
nextflow run ../main.nf -stub-run -params-file stub_inputs/pipeline_input.json
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

`concatenate_cpus` (default `16`), `concatenate_memory_gb` (default `64`), and `dragen_scratch_tb` (default `2`) have defaults in `nextflow.config` and don't need to be passed unless you want to override them.

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

### Command-Line Options

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
- **`supersample_qc/`**: final supersample-level report
  - `<supersample_basename>.supersample_qc_metrics.tsv`
  - `<supersample_basename>.guide_assignment_distribution.png` (only if guide assignment ran)
- **`pipeline_info/`**: Nextflow reports (`timeline.html`, `report.html`, `trace.txt`, `dag.svg`)

## CRISPR Guide Assignment

Guide assignment is **enabled by default** (`--run_guide_assignment true`) and adds one step on top of the concatenation that always happens:

1. **Concatenate**: merge all subsamples and extract CRISPR Direct Capture features into `<supersample_basename>.crispr.h5ad` (`bin/concatenate_samples.py`) — this always runs.
2. **Guide assignment**: run CRISPAT's Poisson-Gaussian mixture model on the CRISPR AnnData (`bin/run_guide_assignment.py`) — only if enabled.
3. **Supersample report**: fold guide assignment results into the final QC report (`bin/generate_supersample_qc.py`) — always runs, with or without guide assignment data.

Disable guide assignment with:
```bash
nextflow run main.nf ... --run_guide_assignment false
```

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

**Error**: Concatenation fails with "We cannot process more than 300 guides"
- **Solution**: `bin/concatenate_samples.py` hard-caps CRISPR feature count at 300 for runtime reasons; reduce your guide library or contact the pipeline maintainer if you need this raised.

**Error**: Pipeline exits immediately with "ERROR: --xyz is required"
- **Solution**: Required params are validated by hand-written checks near the top of `main.nf`'s `workflow {}` block, not just the schema — check the exact list of required flags above.

**Out of memory / timeout**: Adjust resource allocations in `nextflow.config` for specific processes, e.g.:
  ```groovy
  process {
      withName: GUIDE_ASSIGNMENT {
          memory = 64.GB
          time = 24.h
      }
  }
  ```

## Additional Resources

- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
- [Nextflow Container Documentation](https://www.nextflow.io/docs/latest/container.html)
- [CRISPAT Documentation](https://github.com/pinellolab/CRISPAT)
