# MDL QC Nextflow Pipeline

This is a pipeline for processing single-cell QC metrics using Nextflow.

## Overview

This pipeline ingests metrics CSVs and data files, processes them with a Python script, and generates output CSV files.

## Pipeline Structure

```
mdl_qc/qc_pipeline/
├── main.nf                 # Main Nextflow pipeline
├── nextflow.config         # Pipeline configuration
├── modules/
│   └── process.nf          # Process module definitions
├── bin/
│   └── process_metrics.py  # Python processing script (implement your logic here)
└── README.md               # This file
```

## Quick Start

### Prerequisites

- Nextflow (>= 22.10.0)
- Python 3.8+
- pandas

Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
```

### Running the Pipeline

Basic usage:
```bash
nextflow run main.nf \
  --input_metrics 'data/metrics/*.csv' \
  --input_data 'data/additional_data.csv' \
  --outdir results
```

### Command-Line Options

- `--input_metrics`: Path to input metrics CSV file(s). Supports glob patterns (e.g., `'data/*.csv'`)
- `--input_data`: Path to additional input data files
- `--outdir`: Output directory (default: `results`)
- `--help`: Show help message

### Execution Profiles

Run with different executors:

```bash
# Local execution without containers (default)
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --scrna_metrics data/metrics.csv \
  --data_filtered_matrix data/matrix.mtx \
  --data_filtered_barcodes data/barcodes.tsv \
  --data_filtered_features data/features.tsv

# AWS Batch with Docker containers
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --scrna_metrics s3://bucket/metrics.csv \
  --data_filtered_matrix s3://bucket/matrix.mtx \
  --data_filtered_barcodes s3://bucket/barcodes.tsv \
  --data_filtered_features s3://bucket/features.tsv \
  --outdir s3://bucket/results

# SLURM cluster
nextflow run main.nf -profile cluster --num_input_cells ...

# Google Cloud Platform with containers
nextflow run main.nf -profile gcp --num_input_cells ...

# Local with Docker
nextflow run main.nf -profile docker --num_input_cells ...
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

## Implementing Your Logic

The placeholder Python script is located at `bin/process_metrics.py`. Implement your processing logic in these functions:

1. **`load_metrics()`** - Load and parse metrics CSV files
2. **`load_data()`** - Load additional data files
3. **`process_data()`** - Main processing logic (merge, calculate, transform)
4. **`save_output()`** - Save results (already implemented)

The script receives three arguments:
- `--metrics`: Path to metrics CSV
- `--data`: Path to data file
- `--output`: Path for output CSV

## Output

Results are published to `${params.outdir}/processed/`:
- `*.output.csv` - Processed results
- `*.log` - Processing logs

Pipeline execution reports are saved to `${params.outdir}/pipeline_info/`:
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
     params.container_metrics = 'account.dkr.ecr.region.amazonaws.com/image:tag'
     params.container_guide_assignment = 'account.dkr.ecr.region.amazonaws.com/image:tag'
     ```
   - Or via command line:
     ```bash
     --container_metrics 'your-ecr-uri:tag' \
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
   pip install pandas numpy scipy  # Add other required packages
   ```

2. **Ensure scripts are executable**:
   ```bash
   chmod +x bin/process_metrics.py
   chmod +x bin/run_guide_assignment.py
   ```

3. **Run with local profile**:
   ```bash
   nextflow run main.nf -profile local --num_input_cells ...
   ```

### Additional Configuration

Modify `nextflow.config` to adjust:
- Resource allocations (CPU, memory, time)
- Executor settings
- Process-specific configurations
- GCP project settings (for cloud execution)

## Example Workflow

1. Place your metrics CSVs in a directory (e.g., `data/metrics/`)
2. Place your additional data files in a location (e.g., `data/reference.csv`)
3. Implement your processing logic in `bin/process_metrics.py`
4. Run the pipeline:
   ```bash
   nextflow run main.nf \
     --input_metrics 'data/metrics/*.csv' \
     --input_data 'data/reference.csv' \
     --outdir results
   ```
5. Check outputs in `results/processed/`

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run main.nf --input_metrics ... -resume
```

## Troubleshooting

- **Permission denied**: Make sure `bin/process_metrics.py` is executable: `chmod +x bin/process_metrics.py`
- **Module not found**: Install required Python packages in your environment
- **Out of memory**: Adjust memory allocation in `nextflow.config` under `process.withName.PROCESS_METRICS`