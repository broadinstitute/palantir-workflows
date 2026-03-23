# Dockerfile for PIPseq QC Pipeline
FROM continuumio/miniconda3:24.1.2-0

# Set working directory
WORKDIR /app

# Install system dependencies and AWS CLI
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI v2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf aws awscliv2.zip

# Create conda environment with Python packages
RUN conda install -y -c conda-forge -c bioconda \
    python=3.10 \
    pandas \
    numpy \
    scanpy \
    anndata \
    matplotlib \
    scipy \
    && conda clean -a -y

    # Copy requirements.txt and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install Node.js 24 from official prebuilt binaries and puppeteer
RUN curl -fsSL https://nodejs.org/dist/v24.14.0/node-v24.14.0-linux-x64.tar.xz | tar -xJ -C /usr/local --strip-components=1 && \
    npm install -g puppeteer

RUN git clone https://github.com/michaelgatzen/sankeymatic_local.git

# Set PATH for scripts
ENV PATH="/app/bin:${PATH}"

CMD ["/bin/bash"]
