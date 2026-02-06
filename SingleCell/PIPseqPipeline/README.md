# PIPseq QC Nextflow Pipeline

This is a pipeline for processing single-cell QC metrics from PIPseq data using Nextflow.

## Overview

This pipeline processes single-cell RNA-seq data with optional CRISPR guide assignment. It extracts CRISPR features from subsamples (which can be concatenated into a supersample), performs guide assignment using CRISPAT, and generates comprehensive QC reports.

### Terminology

- **Subsample**: An individual sample (e.g., a technical replicate or a well from a plate)
- **Supersample**: A logical grouping of subsamples that should be analyzed together (e.g., all wells from one biological sample)

### Pipeline Workflow

The pipeline performs three main tasks:
1. **Concatenate Subsamples**: Combines subsamples into a supersample (handles single subsample case automatically)
2. **Extract CRISPR Features**: Extracts "CRISPR Direct Capture" features from the concatenated data
3. **Guide Assignment**: Performs CRISPAT guide assignment using Poisson-Gaussian mixture model (optional)
4. **Generate Report Data**: Creates comprehensive QC metrics and visualizations for each subsample

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
└── README.md                        # This file
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

```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --samplesheet samplesheet.csv \
  --supersample_id "Supersample_A" \
  --supersample_basename "supersample_A" \
  --outdir results
```

Example samplesheet.csv:
```csv
data_dir,dragen_file_prefix,subsample_id,subsample_basename
/path/to/data1,sample_001,Subsample_001,subsample_001
/path/to/data2,sample_002,Subsample_002,subsample_002
/path/to/data3,sample_003,Subsample_003,subsample_003
```

**Pipeline behavior:**
- Concatenates all subsamples (handles single subsample case automatically)
- Runs GUIDE_ASSIGNMENT on concatenated CRISPR features (if enabled)
- Per-subsample QC reports are always generated in `outdir/<supersample_basename>/<subsample_basename>/qc/`
- Concatenated AnnData outputs to `outdir/<supersample_basename>/adata/`
- Guide assignments are output to `outdir/<supersample_basename>/crispat_ga/`

### Command-Line Options

**Required:**
- `--num_input_cells`: Number of input cells (integer)
- `--samplesheet`: CSV file with columns: data_dir, dragen_file_prefix, subsample_id, subsample_basename
- `--supersample_id`: Supersample identifier
- `--supersample_basename`: Supersample basename for output organization

**Optional:**
- `--run_guide_assignment`: Run CRISPR guide assignment (default: true)
- `--outdir`: Output directory (default: `results`)
- `--help`: Show help message

**Expected files in each data_dir:**
- `<dragen_file_prefix>.scRNA_metrics.csv`
- `<dragen_file_prefix>.scRNA.barcodeSummary.tsv`
- `<dragen_file_prefix>.scRNA.filtered.matrix.mtx.gz`
- `<dragen_file_prefix>.scRNA.filtered.barcodes.tsv.gz`
- `<dragen_file_prefix>.scRNA.filtered.features.tsv.gz`

### Execution Profiles

Run with different executors:

```bash
# Local execution without containers (default)
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --samplesheet samplesheet.csv \
  --outdir results

# AWS Batch with Docker containers
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --samplesheet s3://bucket/samplesheet.csv \
  --outdir s3://bucket/results

# SLURM cluster
nextflow run main.nf -profile cluster \
  --num_input_cells 10000 \
  --samplesheet samplesheet.csv \
  --outdir results

# Google Cloud Platform with containers
nextflow run main.nf -profile gcp \
  --num_input_cells 10000 \
  --samplesheet gs://bucket/samplesheet.csv \
  --outdir gs://bucket/results

# Local with Docker
nextflow run main.nf -profile docker \
  --num_input_cells 10000 \
  --samplesheet samplesheet.csv \
  --outdir results
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

**`${params.outdir}/<supersample_basename>/<subsample_basename>/qc/`:** (per-subsample QC reports)
- `<subsample_basename>.qc_barcode_metrics.tsv` - Barcode rank metrics
- `<subsample_basename>.qc_metrics.tsv` - Single-cell RNA QC metrics
- `<subsample_basename>.qc_guide_assignment_stats.tsv` - Guide assignment statistics (if enabled)
- `<subsample_basename>.qc_guide_assignment_distribution.png` - Guide assignment distribution plot (if enabled)

**`${params.outdir}/<supersample_basename>/adata/`:** (concatenated AnnData files)
- `<supersample_basename>.h5ad` - Concatenated AnnData with all features
- `<supersample_basename>.crispr.h5ad` - Concatenated CRISPR features only

**`${params.outdir}/<supersample_basename>/crispat_ga/`:** (guide assignments - if enabled)
- `assignments.csv` - CRISPAT guide assignments

**`${params.outdir}/pipeline_info/`:** (Nextflow reports)
- `timeline.html` - Execution timeline
- `report.html` - Resource usage report
- `trace.txt` - Detailed execution trace
- `dag.svg` - Pipeline DAG visualization

## Configuration

### Resource Allocation

Modify `nextflow.config` to adjust:
- Resource allocations (CPU, memory, time)
- Executor settings
- Process-specific configurations (EXTRACT_CRISPR_FEATURES, GUIDE_ASSIGNMENT, GENERATE_REPORT_DATA)
- Cloud execution settings (AWS Batch, GCP)

Example process-specific configuration:
```groovy
process {
    withName: EXTRACT_CRISPR_FEATURES {
        cpus = 8
        memory = 32.GB
        time = 8.h
    }
    
    withName: GUIDE_ASSIGNMENT {
        cpus = 16
        memory = 64.GB
        time = 24.h
    }
}

### Local Execution Setup

**Use when**: Development, testing, or when you have Python dependencies installed locally.

**Requirements**:
- Python 3.8+ with pandas, scanpy, anndata, crispat, matplotlib, numpy, scipy installed
- Nextflow >= 22.10.0
- Scripts in `bin/` directory must be executable

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
     --dragen_file_prefix sample_123 \
     --sample_id "Sample_123" \
     --sample_basename "sample_123" \
     --outdir results
   ```

### AWS Batch Execution Setup

**Use when**: Production runs, large-scale processing, cloud-native workflows.

**Prerequisites**:
1. **AWS Batch Setup**: Create compute environment and job queue
2. **ECR Repositories**: Push Docker images to ECR
3. **S3 Buckets**: Create buckets for data and work directory
4. **IAM Permissions**: Ensure proper S3 access for Batch execution role

**Configuration**:

Update `nextflow.config` with your AWS details:
```groovy
awsbatch {
    aws.region = 'us-east-1'                    // Your region
    workDir = 's3://my-bucket/work'             // S3 work directory
    process.queue = 'my-batch-queue'            // Batch queue name
    
    // Container images in ECR
    params.container_metrics = '123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-metrics:v1.0'
    params.container_crispr = '123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-crispr:v1.0'
    params.container_guide_assignment = '123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-guide:v1.0'
}
```

**Execution**:
```bash
# Configure AWS credentials
export AWS_ACCESS_KEY_ID=your_key_id
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_DEFAULT_REGION=us-east-1

# Run pipeline with S3 paths
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --data_dir s3://my-bucket/data/sample_output \
  --dragen_file_prefix sample_123 \
  --sample_id "Sample_123" \
  --sample_basename "sample_123" \
  --outdir s3://my-bucket/results
```

**AWS Batch Tips**:
- **Resume failed runs**: Add `-resume` flag
- **Monitor**: Check AWS Batch console for job status
- **Logs**: CloudWatch Logs for each task
- **Costs**: S3 storage + EC2 compute time
- **Data transfer**: Input/output to S3 is fastest

### Additional Configuration

Modify `nextflow.config` to adjust:
- Resource allocations (CPU, memory, time)
- Executor settings
- Process-specific configurations (EXTRACT_CRISPR_FEATURES, GUIDE_ASSIGNMENT, GENERATE_REPORT_DATA)
- Cloud execution settings (AWS Batch, GCP)

## Example Workflow

1. Prepare your samplesheet with the required data:
   ```csv
   data_dir,dragen_file_prefix,sample_id,sample_basename
   /path/to/data,sample_123,Sample_123,sample_123
   ```
   
   Each data_dir should contain:
   ```
   data/sample_output/
   ├── sample_123.scRNA_metrics.csv
   ├── sample_123.scRNA.barcodeSummary.tsv
   ├── sample_123.scRNA.filtered.matrix.mtx.gz
   ├── sample_123.scRNA.filtered.barcodes.tsv.gz
   └── sample_123.scRNA.filtered.features.tsv.gz
   ```

2. Run the pipeline:
   ```bash
   nextflow run main.nf \
     --num_input_cells 10000 \
     --samplesheet samplesheet.csv \
     --outdir results
   ```

3. Check outputs in the results directory:
   - QC reports per sample: `results/<sample_basename>/qc/`
   - Guide assignments: `results/ga_crispat/`
   - CRISPR h5ad (single sample): `results/crispr_adata/`
   - Concatenated h5ad (multiple samples): `results/concatenated/`

## CRISPR Guide Assignment

Guide assignment is **enabled by default** (`--run_guide_assignment true`).

**Three-step process (single sample)**:
1. **Extract CRISPR features**: Reads matrix/barcodes/features files, extracts "CRISPR Direct Capture" feature types, and writes to an h5ad file (`bin/extract_crispr_features.py`)
2. **Run guide assignment**: Uses the h5ad file to run CRISPAT's Poisson-Gaussian mixture model (`bin/run_guide_assignment.py`)
3. **Generate report**: Combines all data to create comprehensive QC metrics (`bin/generate_report_data.py`)

**Multi-sample process**:
1. **Concatenate samples**: Combines all samples into a single AnnData object (`bin/concatenate_samples.py`)
2. **Run guide assignment**: Uses the concatenated h5ad to run CRISPAT
3. **Generate reports**: Creates per-sample QC metrics

**Disable guide assignment**:
```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --samplesheet samplesheet.csv \
  --run_guide_assignment false \
  --outdir results
```

**Output files**:
- `crispr_adata/<sample_basename>.crispr.h5ad`: Extracted CRISPR features (single sample only)
- `concatenated/concatenated.h5ad`: Concatenated CRISPR features (multiple samples)
- `ga_crispat/<sample_basename>.assignments.csv` or `ga_crispat/concatenated.assignments.csv`: Guide assignments
- `<sample_basename>/qc/<sample_basename>.qc_*.tsv`: Per-sample QC metrics with guide assignment data integrated

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run main.nf \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --dragen_file_prefix sample_123 \
  --sample_id "Sample_123" \
  --sample_basename "sample_123" \
  -resume
```

## Profile Comparison

| Profile | Executor | Containers | Best For |
|---------|----------|------------|----------|
| `local` | Local machine | None | Development, testing, small datasets |
| `awsbatch` | AWS Batch | Docker (ECR) | Production, large-scale, cloud-native |
| `docker` | Local machine | Docker | Local testing with containers |
| `cluster` | SLURM | Optional | HPC clusters |
| `gcp` | Google Cloud | Docker | Google Cloud users |

## Troubleshooting

### Local Profile Issues

**Error**: `command not found: extract_crispr_features.py`
- **Solution**: Ensure scripts are executable: `chmod +x bin/*.py`
- **Solution**: Check Python is in PATH

**Error**: `ModuleNotFoundError: No module named 'scanpy'` (or anndata, crispat)
- **Solution**: Install dependencies: `pip install pandas scanpy anndata crispat matplotlib numpy scipy`

**Error**: CRISPR features extraction fails
- **Solution**: Verify your features file contains "CRISPR Direct Capture" feature type entries

**Error**: File not found
- **Solution**: Verify file paths and ensure the `dragen_file_prefix` matches your input files exactly

### AWS Batch Profile Issues

**Error**: `Unable to locate credentials`
- **Solution**: Configure AWS CLI or set environment variables

**Error**: `Access Denied` for S3
- **Solution**: Check IAM role attached to AWS Batch compute environment

**Error**: Container fails to start
- **Solution**: Verify ECR image URI is correct
- **Solution**: Check Batch execution role has ECR pull permissions

**Error**: Job stays in RUNNABLE state
- **Solution**: Check compute environment is ENABLED and VALID
- **Solution**: Verify VPC/subnet configuration

### General Issues

- **Out of memory**: Adjust memory allocation in `nextflow.config` for specific processes:
  ```groovy
  withName: EXTRACT_CRISPR_FEATURES {
      memory = 32.GB
  }
  ```
- **Guide assignment fails**: Check that CRISPR Direct Capture features are present in your features file

## Additional Resources

- [Nextflow AWS Batch Documentation](https://www.nextflow.io/docs/latest/awscloud.html)
- [AWS Batch Setup Guide](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html)
- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
- [CRISPAT Documentation](https://github.com/pinellolab/CRISPAT)
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