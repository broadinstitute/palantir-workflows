# This Dockerfile is used to build the monolithic docker image used for most of the tasks in the BenchmarkSVs and CleanSVs workflows.
# It contains htslib, bcftools, truvari, pandas, bedtools, mafft, hiphase.
FROM ubuntu:24.04
RUN apt-get update && apt-get install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev
# RUN pip install --break-system-packages --upgrade pip
RUN pip install --break-system-packages pywfa==0.5.1 truvari==4.3.1 pandas==2.2.3 pysam==0.22.1

# Install htslib
RUN git clone --recurse-submodules --branch 1.21 https://github.com/samtools/htslib.git
WORKDIR htslib
RUN make && make install

# Install bcftools
WORKDIR /
RUN git clone --branch 1.21 https://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
ENV BCFTOOLS_PLUGINS="/bcftools/plugins/"

## Add bedtools
WORKDIR /
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static .
RUN mv bedtools.static bedtools && chmod +x bedtools && cp bedtools /bin/

## Add MAFFT (for truvari refine step)
WORKDIR /
RUN git clone https://github.com/GSLBiotech/mafft
WORKDIR mafft/core
RUN make install
WORKDIR /

## Add HiPhase (to phase SVs before truvari refine)
RUN mkdir -p hiphase
ADD https://github.com/PacificBiosciences/HiPhase/releases/download/v1.4.5/hiphase-v1.4.5-x86_64-unknown-linux-gnu.tar.gz hiphase
WORKDIR hiphase
RUN tar -xvzf hiphase-v1.4.5-x86_64-unknown-linux-gnu.tar.gz
RUN cp hiphase-v1.4.5-x86_64-unknown-linux-gnu/hiphase /bin
WORKDIR /
