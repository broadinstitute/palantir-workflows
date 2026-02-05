# Dockerfile for CRISPR Guide Assignment
FROM continuumio/miniconda3:24.1.2-0

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY docker/guide_assignment/requirements.txt .

# Create conda environment with Python packages for guide assignment
RUN conda install -y -c conda-forge -c bioconda \
    python=3.10 \
    pandas \
    numpy \
    scanpy \
    anndata \
    crispat \
    pyro-ppl \
    && conda clean -a -y

# Install additional pip packages if file exists
RUN if [ -s requirements.txt ]; then pip install --no-cache-dir -r requirements.txt; fi

# Make scripts executable
RUN mkdir -p /app/bin

# Set PATH
ENV PATH="/app/bin:${PATH}"

CMD ["/bin/bash"]
