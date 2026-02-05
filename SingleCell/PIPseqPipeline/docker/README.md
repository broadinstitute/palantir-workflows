# Docker and Container Management for Nextflow Pipeline

## Overview

Similar to WDL's `runtime.docker` attribute, Nextflow allows you to specify Docker containers at the **process level** or **globally** in the config file.

This pipeline uses **separate Docker images** for different processes:
- **Metrics Processing**: Lighter image with core dependencies
- **Guide Assignment**: Specialized image with CRISPAT and related tools

## Directory Structure

```
docker/
├── metrics/
│   ├── Dockerfile           # Metrics processing image
│   └── requirements.txt     # Python pip dependencies
└── guide_assignment/
    ├── Dockerfile           # CRISPR guide assignment image
    └── requirements.txt     # Python pip dependencies
```

## Quick Comparison: WDL vs Nextflow

**WDL:**
```wdl
task my_task {
    runtime {
        docker: "my-image:latest"
    }
}
```

**Nextflow:**
```groovy
process MY_PROCESS {
    container "my-image:latest"
    
    script:
    """
    # your code
    """
}
```

## Setup for Illumina Connected Analytics (ICA)

### 1. Build and Push Your Docker Images

**Metrics Processing Image:**
```bash
cd mdl_qc/qc_pipeline
docker build -f docker/metrics/Dockerfile -t your-dockerhub-username/singlecell-qc-metrics:latest .

# Test locally (optional)
docker run -it your-dockerhub-username/singlecell-qc-metrics:latest python --version

# Push to Docker Hub (or your registry)
docker push your-dockerhub-username/singlecell-qc-metrics:latest
```

**Guide Assignment Image:**
```bash
docker build -f docker/guide_assignment/Dockerfile -t your-dockerhub-username/singlecell-qc-guide-assignment:latest .

# Test locally (optional)
docker run -it your-dockerhub-username/singlecell-qc-guide-assignment:latest python -c "import crispat; print('crispat loaded')"

# Push to Docker Hub (or your registry)
docker push your-dockerhub-username/singlecell-qc-guide-assignment:latest
```

**Build both at once:**
```bash
# Build both images
docker build -f docker/metrics/Dockerfile -t your-dockerhub-username/singlecell-qc-metrics:latest .
docker build -f docker/guide_assignment/Dockerfile -t your-dockerhub-username/singlecell-qc-guide-assignment:latest .

# Push both
docker push your-dockerhub-username/singlecell-qc-metrics:latest
docker push your-dockerhub-username/singlecell-qc-guide-assignment:latest
```

**Alternative: Use AWS ECR for ICA**
```bash
# Login to AWS ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin <account-id>.dkr.ecr.us-east-1.amazonaws.com

# Create repositories
aws ecr create-repository --repository-name singlecell-qc-metrics --region us-east-1
aws ecr create-repository --repository-name singlecell-qc-guide-assignment --region us-east-1

# Tag for ECR
docker tag singlecell-qc-metrics:latest <account-id>.dkr.ecr.us-east-1.amazonaws.com/singlecell-qc-metrics:latest
docker tag singlecell-qc-guide-assignment:latest <account-id>.dkr.ecr.us-east-1.amazonaws.com/singlecell-qc-guide-assignment:latest

# Push to ECR
docker push <account-id>.dkr.ecr.us-east-1.amazonaws.com/singlecell-qc-metrics:latest
docker push <account-id>.dkr.ecr.us-east-1.amazonaws.com/singlecell-qc-guide-assignment:latest
```

### 2. Run on ICA

**Option A: Specify containers in command line**
```bash
nextflow run main.nf \
  -profile ica \
  --container-metrics your-dockerhub-username/singlecell-qc-metrics:latest \
  --container-guide-assignment your-dockerhub-username/singlecell-qc-guide-assignment:latest \
  --num-input-cells 10000 \
  --scrna-metrics s3://bucket/metrics.csv \
  --data-filtered-matrix s3://bucket/matrix.mtx.gz \
  --data-filtered-barcodes s3://bucket/barcodes.tsv.gz \
  --data-filtered-features s3://bucket/features.tsv.gz \
  --feature-barcode-reference s3://bucket/barcode_ref.csv
```

**Option B: Update process definitions**

Edit `modules/process.nf` to replace:
```groovy
container "${params.container_metrics ?: 'your-dockerhub-username/singlecell-qc-metrics:latest'}"
```

Edit `modules/guide_assignment.nf` to replace:
```groovy
container "${params.container_guide_assignment ?: 'your-dockerhub-username/singlecell-qc-guide-assignment:latest'}"
```
with your actual image names.

### 3. ICA-Specific Configuration

The pipeline includes an `ica` profile in `nextflow.config` with AWS Batch settings. Key points:

- **Docker enabled by default** when using `-profile ica`
- **AWS Batch executor** for cloud execution
- **S3 work directory** for intermediate files
- **Region**: Update `aws.region` in config to match your ICA deployment

## Container Management Strategy (Current Setup)

This pipeline uses **process-specific containers** - the optimal approach for modularity:

```groovy
// In modules/process.nf
process PROCESS_METRICS {
    container "your-username/singlecell-qc-metrics:latest"
    // ...
}

// In modules/guide_assignment.nf
process GUIDE_ASSIGNMENT {
    container "your-username/singlecell-qc-guide-assignment:latest"
    // ...
}
```

### Benefits of Separate Containers

✅ **Smaller images**: Metrics processing doesn't need CRISPAT  
✅ **Faster builds**: Changes to guide assignment don't rebuild metrics image  
✅ **Better isolation**: Different dependency versions if needed  
✅ **Clearer dependencies**: Each process has exactly what it needs  

## Updating the Dockerfiles

### Metrics Processing (`docker/metrics/Dockerfile`)
Includes:
- ✅ Python 3.10
- ✅ pandas, numpy, scanpy, anndata
- ✅ python-slugify
- ✅ Additional packages from `docker/metrics/requirements.txt`

### Guide Assignment (`docker/guide_assignment/Dockerfile`)
Includes:
- ✅ Python 3.10
- ✅ pandas, numpy, scanpy, anndata
- ✅ **crispat** (CRISPR-specific)
- ✅ **pyro-ppl** (probabilistic programming for CRISPAT)
- ✅ Additional packages from `docker/guide_assignment/requirements.txt`

### Adding Packages

**To metrics processing:**
1. **Conda packages** (preferred for scientific packages):
   ```dockerfile
   # Edit docker/metrics/Dockerfile
   RUN conda install -y -c conda-forge -c bioconda \
       your-package-name \
       && conda clean -a -y
   ```

2. **Pip packages**:
   Edit `docker/metrics/requirements.txt` and rebuild

**To guide assignment:**
1. Edit `docker/guide_assignment/Dockerfile` or `docker/guide_assignment/requirements.txt`
2. Rebuild only the guide assignment image

## Testing Locally

Before pushing to ICA, test locally:

```bash
# Run with Docker profile
nextflow run main.nf \
  -profile docker \
  --container-metrics your-dockerhub-username/singlecell-qc-metrics:latest \
  --container-guide-assignment your-dockerhub-username/singlecell-qc-guide-assignment:latest \
  --num-input-cells 1000 \
  --scrna-metrics test/data/metrics.csv \
  --data-filtered-matrix test/data/matrix.mtx.gz \
  --data-filtered-barcodes test/data/barcodes.tsv.gz \
  --data-filtered-features test/data/features.tsv.gz
```

**Test individual containers:**
```bash
# Test metrics container
docker run -it your-dockerhub-username/singlecell-qc-metrics:latest python -c "import pandas, scanpy; print('OK')"

# Test guide assignment container
docker run -it your-dockerhub-username/singlecell-qc-guide-assignment:latest python -c "import crispat; print('OK')"
```

## Troubleshooting

**Problem: "Container not found"**
- Ensure image is pushed to registry accessible by ICA
- Check image name/tag spelling

**Problem: "Permission denied" for scripts**
- Scripts in `bin/` are automatically made available in PATH
- No need to `chmod +x` inside container

**Problem: "Module not found" in Python**
- Add missing package to Dockerfile
- Rebuild and push new image

**Problem: Different behavior locally vs ICA**
- Check platform compatibility: add `--platform linux/amd64` to docker build
- Verify AWS region matches ICA deployment

## Additional Resources

- Nextflow containers: https://www.nextflow.io/docs/latest/container.html
- ICA documentation: Check Illumina's ICA user guide for AWS Batch configuration
- AWS Batch: https://www.nextflow.io/docs/latest/awscloud.html#aws-batch
