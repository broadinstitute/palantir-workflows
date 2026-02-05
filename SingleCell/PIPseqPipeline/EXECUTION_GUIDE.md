# PIPseq Pipeline Execution Guide

Quick reference for running the PIPseq single-cell QC pipeline in different environments.

## Pipeline Overview

The pipeline performs three main tasks:
1. **Extract CRISPR Features**: Reads 10x matrix files and extracts "CRISPR Direct Capture" features to h5ad format
2. **Guide Assignment**: Performs CRISPAT guide assignment using Poisson-Gaussian mixture model (optional)
3. **Generate Report Data**: Creates comprehensive QC metrics and visualizations

## Local Execution (No Containers)

**Use when**: Development, testing, or when you have Python dependencies installed locally.

```bash
# Make sure Python scripts are executable
chmod +x bin/*.py

# Run the pipeline
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  --outdir results
```

**Requirements**:
- Python 3.8+ with pandas, scanpy, anndata, crispat, matplotlib, numpy, scipy installed
- Nextflow >= 22.10.0
- Scripts in `bin/` directory must be executable

**Expected files in data_dir**: 
- `<prefix>.scRNA_metrics.csv`
- `<prefix>.filtered.matrix.mtx.gz`
- `<prefix>.filtered.barcodes.tsv.gz`
- `<prefix>.filtered.features.tsv.gz`

---

## AWS Batch Execution (Docker Containers)

**Use when**: Production runs, large-scale processing, cloud-native workflows.

### Prerequisites

1. **AWS Batch Setup**: Create compute environment and job queue
2. **ECR Repositories**: Push Docker images to ECR
3. **S3 Buckets**: Create buckets for data and work directory
4. **IAM Permissions**: Ensure proper S3 access for Batch execution role

### Configuration

Update [nextflow.config](nextflow.config) with your AWS details:

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

### Execution

```bash
# Configure AWS credentials
export AWS_ACCESS_KEY_ID=your_key_id
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_DEFAULT_REGION=us-east-1

# Run pipeline with S3 paths
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --data_dir s3://my-bucket/data/sample_output \
  --file_prefix sample_123 \
  --outdir s3://my-bucket/results

# Optional: Override container images at runtime
nextflow run main.nf -profile awsbatch \
  --container_metrics '123456789012.dkr.ecr.us-east-1.amazonaws.com/custom:tag' \
  --container_guide_assignment '123456789012.dkr.ecr.us-east-1.amazonaws.com/custom:tag' \
  --scrna_metrics s3://...
```

### AWS Batch Tips

- **Resume failed runs**: Add `-resume` flag
- **Monitor**: Check AWS Batch console for job status
- **Logs**: CloudWatch Logs for each task
- **Costs**: S3 storage + EC2 compute time
- **Data transfer**: Input/output to S3 is fastest

---

## Optional: CRISPR Guide Assignment

Enable guide assignment with `--run_guide_assignment true` (enabled by default).

**Three-step process**:
1. **Extract CRISPR features**: Reads matrix/barcodes/features files, extracts "CRISPR Direct Capture" feature types, and writes to an h5ad file (`bin/extract_crispr_features.py`)
2. **Run guide assignment**: Uses the h5ad file to run CRISPAT's Poisson-Gaussian mixture model (`bin/run_guide_assignment.py`)
3. **Generate report**: Combines all data to create comprehensive QC metrics (`bin/generate_report_data.py`)

```bash
# Local with guide assignment (default)
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  --outdir results

# AWS Batch with guide assignment
nextflow run main.nf -profile awsbatch \
  --num_input_cells 10000 \
  --data_dir s3://bucket/sample_output \
  --file_prefix sample_123 \
  --outdir s3://bucket/results

# Disable guide assignment
nextflow run main.nf -profile local \
  --num_input_cells 10000 \
  --data_dir data/sample_output \
  --file_prefix sample_123 \
  --run_guide_assignment false \
  --outdir results
```

**Output files**:
- `crispr_adata/<prefix>.crispr.h5ad`: Extracted CRISPR features in AnnData format
- `ga_crispat/<prefix>.assignments.csv`: Guide assignments from CRISPAT
- `processed/<prefix>.output.csv`: Final QC metrics with guide assignment data integrated

---

## Profile Comparison

| Profile | Executor | Containers | Best For |
|---------|----------|------------|----------|
| `local` | Local machine | None | Development, testing, small datasets |
| `awsbatch` | AWS Batch | Docker (ECR) | Production, large-scale, cloud-native |
| `docker` | Local machine | Docker | Local testing with containers |
| `cluster` | SLURM | Optional | HPC clusters |
| `gcp` | Google Cloud | Docker | Google Cloud users |

---

## Troubleshooting

### Local Profile Issues

**Error**: `command not found: extract_crispr_features.py`
- **Solution**: Ensure scripts are executable: `chmod +x bin/*.py`
- **Solution**: Check Python is in PATH

**Error**: `ModuleNotFoundError: No module named 'scanpy'` (or anndata, crispat)
- **Solution**: Install dependencies: `pip install pandas scanpy anndata crispat matplotlib numpy scipy`

**Error**: CRISPR features extraction fails
- **Solution**: Verify your features file contains "CRISPR Direct Capture" feature type entries

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

---

## Additional Resources

- [Nextflow AWS Batch Documentation](https://www.nextflow.io/docs/latest/awscloud.html)
- [AWS Batch Setup Guide](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html)
- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
