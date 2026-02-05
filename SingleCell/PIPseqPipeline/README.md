# PIPseq QC Nextflow Pipeline

This is a pipeline for processing single-cell QC metrics from PIPseq data using Nextflow.

## Overview

This pipeline processes single-cell RNA-seq data with optional CRISPR guide assignment. It extracts CRISPR features, performs guide assignment using CRISPAT, and generates comprehensive QC reports.

## Pipeline Structure

```
SingleCell/PIPseqPipeline/
├── main.nf                          # Main Nextflow pipeline
├── nextflow.config                  # Pipeline configuration
├── modules/
│   ├── extract_crispr_features.nf   # Extract CRISPR Direct Capture features
│   ├── guide_assignment.nf          # CRISPAT guide assignment
│   └── generate_report_data.nf      # Generate QC report data
├── bin/
│   ├── extract_crispr_features.py   # Extract CRISPR features to h5ad
│   ├── run_guide_assignment.py      # Run CRISPAT guide assignment
│   └── generate_report_data.py      # Generate final QC metrics
├── README.md                        # This file
└── EXECUTION_GUIDE.md               # Detailed execution guide
```

## Quick Start

### Prerequisites

- Nextflow (>= 22.10.0)
- Python 3.8+
- Required Python packages: pandas, scanpy, anndata, crispat, matplotlib

Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
```

Install Python dependencies:
```bash
pip install pandas scanpy anndata crispat matplotlib
```

### Running the Pipeline

Basic usage:
```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  --outdir results
```

### Command-Line Options

**Required:**
- `--num_input_cells`: Number of input cells (integer)
- `--data_dir`: Directory containing all input data files
- `--file_prefix`: Common prefix for all input files

**Optional:**
- `--run_guide_assignment`: Run CRISPR guide assignment (default: true)
- `--outdir`: Output directory (default: `results`)
- `--help`: Show help message

**Expected files in data_dir:**
- `<prefix>.scRNA_metrics.csv`
- `<prefix>.filtered.matrix.mtx.gz`
- `<prefix>.filtered.barcodes.tsv.gz`
- `<prefix>.filtered.features.tsv.gz`

### Execution Profiles

Run with different executors:

```bash
# Local execution without containers (default)
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  --outdir results

# AWS Batch with Docker containers
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --data_dir s3://bucket/sample_output \
  --file_prefix sample_123 \
  --outdir s3://bucket/results

# SLURM cluster
nextflow run main.nf -profile cluster \
  --num_input_cells 10000 \
  --data_dir /path/to/data \
  --file_prefix sample_123

# Google Cloud Platform with containers
nextflow run main.nf -profile gcp \
  --num_input_cells 10000 \
  --data_dir gs://bucket/data \
  --file_prefix sample_123

# Local with Docker
nextflow run main.nf -profile docker \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123
```

#### Profile Details

**`local`** - Local execution without containers
- Runs directly on your machine
- No Docker or Singularity required
- Python scripts must be available in your environment
- Best for: Development, testing, small datasets

**`awsbatch`** - AWS Batch with Docker containers
- Runs on AWS Batch compute environment
- Requires Docker containers in ECR or Docker Hub
- Data can be in S3 buckets
- Update configuration in `nextflow.config`:
  - `aws.region`: Your AWS region
  - `workDir`: S3 bucket for intermediate files
  - `process.queue`: Your AWS Batch queue name
  - `params.container_metrics`: ECR image URI
  - `params.container_guide_assignment`: ECR image URI
- Best for: Large-scale production runs, cloud-native workflows

## Pipeline Components

### 1. Extract CRISPR Features (`extract_crispr_features.py`)

**Purpose:** Reads 10x-style filtered matrix files and extracts only the "CRISPR Direct Capture" feature types.

**Inputs:**
- Matrix file (`.mtx.gz`)
- Barcodes file (`.tsv.gz`)
- Features file (`.tsv.gz`)

**Output:**
- CRISPR features h5ad file (AnnData format)

### 2. Guide Assignment (`run_guide_assignment.py`)

**Purpose:** Performs CRISPR guide assignment using CRISPAT's Poisson-Gaussian mixture model.

**Input:**
- CRISPR h5ad file from step 1

**Output:**
- Guide assignments CSV file

### 3. Generate Report Data (`generate_report_data.py`)

**Purpose:** Generates comprehensive QC metrics and report data.

**Inputs:**
- scRNA metrics CSV
- Barcode summary
- Sample ID
- Guide assignments (optional)

**Output:**
- QC report CSV with metrics and visualizations

## Output

Results are organized in the output directory:

**`${params.outdir}/crispr_adata/`:**
- `<prefix>.crispr.h5ad` - CRISPR features in AnnData format
- `<prefix>.extract_crispr_features.log` - Processing log

**`${params.outdir}/ga_crispat/`:** (if guide assignment enabled)
- `<prefix>.assignments.csv` - CRISPAT guide assignments
- `<prefix>.guide_assignment.log` - Processing log

**`${params.outdir}/processed/`:**
- `<prefix>.output.csv` - Final QC metrics and report data
- `<prefix>.generate_report_data.log` - Processing log

**`${params.outdir}/pipeline_info/`:** (Nextflow reports)
- `timeline.html` - Execution timeline
- `report.html` - Resource usage report
- `trace.txt` - Detailed execution trace
- `dag.svg` - Pipeline DAG visualization

## Configuration

### AWS Batch Setup

To use the `awsbatch` profile, you need to:

1. **Update `nextflow.config`** with your AWS settings:
   ```groovy
   aws.region = 'us-east-1'              // Your AWS region
   workDir = 's3://your-bucket/work'     // S3 bucket for work directory
   process.queue = 'your-batch-queue'    // AWS Batch queue name
   ```

2. **Set container images** (choose one option):
   - In `nextflow.config`:
     ```groovy
     params.container_metrics = 'account.dkr.ecr.region.amazonaws.com/pipseq-metrics:tag'
     params.container_crispr = 'account.dkr.ecr.region.amazonaws.com/pipseq-crispr:tag'
     params.container_guide_assignment = 'account.dkr.ecr.region.amazonaws.com/pipseq-guide:tag'
     ```
   - Or via command line:
     ```bash
     --container_metrics 'your-ecr-uri:tag' \
     --container_crispr 'your-ecr-uri:tag' \
     --container_guide_assignment 'your-ecr-uri:tag'
     ```

3. **AWS Credentials**: Ensure AWS credentials are configured
   ```bash
   aws configure
   # or set environment variables:
   export AWS_ACCESS_KEY_ID=your_key
   export AWS_SECRET_ACCESS_KEY=your_secret
   ```

4. **S3 Permissions**: Ensure your AWS Batch execution role has access to:
   - Input data buckets (read)
   - Output bucket (read/write)
   - Work directory bucket (read/write)

### Local Execution Setup

For local execution without containers:

1. **Install Python dependencies**:
   ```bash
   pip install pandas scanpy anndata crispat matplotlib numpy scipy
   ```

2. **Ensure scripts are executable**:
   ```bash
   chmod +x bin/extract_crispr_features.py
   chmod +x bin/run_guide_assignment.py
   chmod +x bin/generate_report_data.py
   ```

3. **Run with local profile**:
   ```bash
   nextflow run main.nf -profile local \
     --num_input_cells 10000 \
     --data_dir data/sample_output \
     --file_prefix sample_123 \
     --outdir results
   ```

### Additional Configuration

Modify `nextflow.config` to adjust:
- Resource allocations (CPU, memory, time)
- Executor settings
- Process-specific configurations (EXTRACT_CRISPR_FEATURES, GUIDE_ASSIGNMENT, GENERATE_REPORT_DATA)
- Cloud execution settings (AWS Batch, GCP)

## Example Workflow

1. Prepare your data directory with the required files:
   ```
   data/sample_output/
   ├── sample_123.scRNA_metrics.csv
   ├── sample_123.filtered.matrix.mtx.gz
   ├── sample_123.filtered.barcodes.tsv.gz
   └── sample_123.filtered.features.tsv.gz
   ```

2. Run the pipeline:
   ```bash
   nextflow run main.nf \
     --num_input_cells 10000 \
     --data_dir data/sample_output \
     --file_prefix sample_123 \
     --outdir results
   ```

3. Check outputs in the results directory:
   - `results/crispr_adata/` - CRISPR features h5ad
   - `results/ga_crispat/` - Guide assignments (if enabled)
   - `results/processed/` - Final QC metrics

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  -resume
```

## Troubleshooting

- **Permission denied**: Make sure Python scripts are executable:
  ```bash
  chmod +x bin/*.py
  ```
- **ModuleNotFoundError**: Install missing Python packages:
  ```bash
  pip install pandas scanpy anndata crispat matplotlib
  ```
- **Out of memory**: Adjust memory allocation in `nextflow.config` for specific processes:
  ```groovy
  withName: EXTRACT_CRISPR_FEATURES {
      memory = 32.GB
  }
  ```
- **File not found**: Verify file paths and ensure the `file_prefix` matches your input files exactly
- **Guide assignment fails**: Check that CRISPR Direct Capture features are present in your features file

For more detailed execution instructions, see [EXECUTION_GUIDE.md](EXECUTION_GUIDE.md).
- **Module not found**: Install required Python packages in your environment
- **Out of memory**: Adjust memory allocation in `nextflow.config` under `process.withName.PROCESS_METRICS`