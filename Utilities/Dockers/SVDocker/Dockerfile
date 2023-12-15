# This Dockerfile is used to build the monolithic docker image used for most of the tasks in the BenchmarkSVs and CleanSVs workflows.
# It contains htslib, bcftools, truvari, pandas, and bedtools.
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev
RUN pip install truvari==4.0.0 pandas==2.1.3 pysam==0.22.0

# Install htslib
RUN git clone --recurse-submodules --branch 1.18 https://github.com/samtools/htslib.git
WORKDIR htslib
RUN make && make install

# Install bcftools
WORKDIR /
RUN git clone --branch 1.18 https://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
ENV BCFTOOLS_PLUGINS="/bcftools/plugins/"

## Add bedtools
WORKDIR /
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static .
RUN mv bedtools.static bedtools && chmod +x bedtools && cp bedtools /bin/
