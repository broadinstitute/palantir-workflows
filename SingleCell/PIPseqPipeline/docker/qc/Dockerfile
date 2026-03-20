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

# Install Node.js globally (not in /root) and puppeteer
ENV NVM_DIR=/usr/local/nvm

RUN mkdir -p $NVM_DIR && \
    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.4/install.sh | bash && \
    . $NVM_DIR/nvm.sh && \
    nvm install 24 && \
    nvm alias default 24 && \
    nvm use default && \
    npm install -g puppeteer && \
    chmod -R a+rX $NVM_DIR && \
    NODE_ACTUAL_VERSION=$(nvm version default) && \
    ln -sf $NVM_DIR/versions/node/$NODE_ACTUAL_VERSION/bin/node /usr/local/bin/node && \
    ln -sf $NVM_DIR/versions/node/$NODE_ACTUAL_VERSION/bin/npm /usr/local/bin/npm && \
    ln -sf $NVM_DIR/versions/node/$NODE_ACTUAL_VERSION/bin/npx /usr/local/bin/npx && \
    # Create a global node_modules directory that all users can access
    mkdir -p /usr/local/lib/node_modules && \
    chmod -R a+rX /usr/local/lib/node_modules && \
    # Link puppeteer to the global location
    ln -sf $NVM_DIR/versions/node/$NODE_ACTUAL_VERSION/lib/node_modules/puppeteer /usr/local/lib/node_modules/puppeteer

# Set NODE_PATH to the accessible global node_modules
ENV NODE_PATH=/usr/local/lib/node_modules

RUN git clone https://github.com/michaelgatzen/sankeymatic_local.git

# Set PATH for scripts
ENV PATH="/app/bin:${PATH}"

CMD ["/bin/bash"]
