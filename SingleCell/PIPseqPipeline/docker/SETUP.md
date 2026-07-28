# Docker Setup

This directory contains the Docker image definition for the QC-side processes of the PIPseq pipeline (`GENERATE_SUBSAMPLE_QC`, `CONCATENATE`, `CRISPAT_GUIDE_ASSIGNMENT`, `PURITY_BASED_GUIDE_ASSIGNMENT`, `GENERATE_SUPERSAMPLE_QC`).

DRAGEN (`DRAGEN_SCRNA`) uses a separate, Illumina-provided container passed via `--dragen_container`; it is not built from anything in this directory, and it runs on Illumina Connected Analytics (ICA), not AWS Batch.

## Structure

```
docker/
├── README.md              # General Docker documentation
├── SETUP.md               # This file - build/push guide
├── build_and_push.sh      # Script to build and push the qc image to ECR
└── qc/
    ├── Dockerfile         # Single unified image for all non-DRAGEN pipeline processes
    └── requirements.txt   # Python package requirements (optional)
```

## Building and Pushing to ECR

### Prerequisites

1. **AWS CLI configured** with appropriate credentials:
   ```bash
   aws configure
   ```

2. **Docker installed** and running

3. **ECR Registry URL** - Find your ECR registry URL in the AWS Console or run:
   ```bash
   aws ecr describe-repositories --region us-east-1
   ```
   Format: `<account-id>.dkr.ecr.<region>.amazonaws.com`

### Build and Push

Set your ECR registry and build:

```bash
# Set your ECR registry URL
export ECR_REGISTRY=123456789012.dkr.ecr.us-east-1.amazonaws.com
export AWS_REGION=us-east-1

# Build and push with default 'latest' tag
./build_and_push.sh

# Or build with a specific tag (e.g., version or commit hash)
./build_and_push.sh v1.0.0
./build_and_push.sh $(git rev-parse --short HEAD)
```

The script will:
1. Authenticate Docker to your ECR registry
2. Create the ECR repository if it doesn't exist
3. Build the Docker image from `qc/Dockerfile`
4. Push the image to ECR

### Use the pushed image

`qc_container` is a required pipeline parameter (there is no config-file default) — pass it explicitly on the command line, or set `params.qc_container` in a config file included via `-c`:

```bash
nextflow run main.nf \
  --qc_container 123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-qc:latest \
  --dragen_container <your DRAGEN container> \
  --num_input_cells 10000 \
  --fastq_list fastq_list.csv \
  --supersample_id "Sample_A" \
  --supersample_basename "sample_a" \
  --outdir results \
  ...
```

## QC Image Contents

The `qc` image includes:
- **Python 3.10** with conda
- **AWS CLI v2** - Required for S3 file operations
- **Python packages**: pandas, numpy, scanpy, anndata, matplotlib, scipy (via conda); crispat and its dependencies (pyro-ppl, torch, etc.) via `pip install` from a git clone
- **Node.js + Puppeteer + sankeymatic_local** - for the (currently unused/disabled) Sankey plot feature in `bin/generate_supersample_qc.py`
- **System tools**: git, build-essential

This single image is used by all non-DRAGEN pipeline processes:
- `GENERATE_SUBSAMPLE_QC`
- `CONCATENATE`
- `CRISPAT_GUIDE_ASSIGNMENT`
- `PURITY_BASED_GUIDE_ASSIGNMENT`
- `GENERATE_SUPERSAMPLE_QC`

## Troubleshooting

**Authentication errors**: Ensure your AWS credentials are configured and have ECR push permissions

**Repository not found**: The script creates the repository automatically, but ensure your IAM role has `ecr:CreateRepository` permission

**Build failures**: Check that the Dockerfile path is correct and all dependencies are available

**Container fails to start / wrong container used**: Double check the `--qc_container`/`--dragen_container` values passed on the command line (or set in your config file) — the pipeline has no default container image for either.
