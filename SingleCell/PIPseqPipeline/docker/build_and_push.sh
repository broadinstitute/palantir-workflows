#!/bin/bash

# Build and push Docker image for PIPseq QC Pipeline to ECR
# Usage: ./build_and_push.sh [tag]
#   tag: Optional tag for the image (default: latest)

set -e

# Configuration
AWS_REGION="${AWS_REGION:-us-east-1}"
ECR_REGISTRY="${ECR_REGISTRY:-}"  # Set this to your ECR registry URL (e.g., 123456789012.dkr.ecr.us-east-1.amazonaws.com)
REPOSITORY_NAME="pipseq-qc"
IMAGE_TAG="${1:-latest}"

# Check if ECR_REGISTRY is set
if [ -z "$ECR_REGISTRY" ]; then
    echo "ERROR: ECR_REGISTRY environment variable is not set"
    echo "Set it to your ECR registry URL, e.g.:"
    echo "  export ECR_REGISTRY=123456789012.dkr.ecr.us-east-1.amazonaws.com"
    exit 1
fi

# Full image name
IMAGE_NAME="${ECR_REGISTRY}/${REPOSITORY_NAME}:${IMAGE_TAG}"

echo "Building Docker image: ${IMAGE_NAME}"
echo "Region: ${AWS_REGION}"
echo ""

# Navigate to docker directory
cd "$(dirname "$0")"

# Authenticate Docker to ECR
echo "Authenticating to ECR..."
aws ecr get-login-password --region ${AWS_REGION} | docker login --username AWS --password-stdin ${ECR_REGISTRY}

# Create ECR repository if it doesn't exist
echo "Ensuring ECR repository exists..."
aws ecr describe-repositories --repository-names ${REPOSITORY_NAME} --region ${AWS_REGION} 2>/dev/null || \
    aws ecr create-repository --repository-name ${REPOSITORY_NAME} --region ${AWS_REGION}

# Build the Docker image
echo "Building Docker image from qc/Dockerfile..."
docker build -t ${IMAGE_NAME} -f qc/Dockerfile ../

# Push to ECR
echo "Pushing image to ECR..."
docker push ${IMAGE_NAME}

echo ""
echo "✓ Successfully built and pushed: ${IMAGE_NAME}"
echo ""
echo "To use this image in your pipeline, update nextflow.config:"
echo "  params.container_qc = '${IMAGE_NAME}'"
