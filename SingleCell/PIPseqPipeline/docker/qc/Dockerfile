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
    chromium \
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
# Set puppeteer cache to system-wide location accessible to all users
ENV PUPPETEER_CACHE_DIR="/usr/local/share/puppeteer"
ENV NODE_PATH="/usr/local/lib/node_modules"

RUN curl -fsSL https://nodejs.org/dist/v24.14.0/node-v24.14.0-linux-x64.tar.xz | tar -xJ -C /usr/local --strip-components=1 && \
    npm install -g puppeteer

RUN git clone https://github.com/michaelgatzen/sankeymatic_local.git

RUN git clone https://github.com/michaelgatzen/crispat.git
RUN cd crispat && pip install .

# Set PATH for scripts
ENV PATH="/app/bin:${PATH}"

# Give all users full access to /app
RUN chmod -R 777 /app

CMD ["/bin/bash"]
