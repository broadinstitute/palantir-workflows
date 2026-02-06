# Docker Setup for AWS Batch

This directory contains Docker image definitions for running the PIPseq QC Pipeline on AWS Batch.

## Structure

```
docker/
├── README.md              # General Docker documentation
├── SETUP.md               # This file - AWS Batch setup guide
├── build_and_push.sh      # Script to build and push images to ECR
└── qc/
    ├── Dockerfile         # Single unified image for all pipeline processes
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

### Update Nextflow Config

After pushing the image, update `nextflow.config`:

```groovy
profiles {
    awsbatch {
        // Update these values
        aws.region = 'us-east-1'
        workDir = 's3://your-bucket-name/work'
        process.queue = 'your-batch-queue-name'
        params.container_qc = '123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-qc:latest'
    }
}
```

## QC Image Contents

The `qc` image includes:
- **Python 3.10** with conda
- **AWS CLI v2** - Required for S3 file operations in AWS Batch
- **Python packages**: pandas, numpy, scanpy, anndata, crispat, matplotlib, scipy
- **System tools**: git, build-essential

This single image is used by all pipeline processes:
- `GENERATE_REPORT_DATA`
- `CONCATENATE`
- `GUIDE_ASSIGNMENT`

## Running on AWS Batch

Once the image is pushed and configured:

```bash
nextflow run main.nf \
  -profile awsbatch \
  --num_input_cells 10000 \
  --samplesheet s3://your-bucket/samplesheet.csv \
  --supersample_id "Sample_A" \
  --supersample_basename "sample_a" \
  --outdir s3://your-bucket/results
```

## Troubleshooting

**Authentication errors**: Ensure your AWS credentials are configured and have ECR push permissions

**Repository not found**: The script creates the repository automatically, but ensure your IAM role has `ecr:CreateRepository` permission

**Build failures**: Check that the Dockerfile path is correct and all dependencies are available

**AWS Batch issues**: Verify:
- Container image URI is correct in `nextflow.config`
- AWS Batch compute environment has internet access (for pulling from ECR)
- IAM roles have appropriate S3 and ECR permissions
- Work directory S3 bucket exists and is accessible
