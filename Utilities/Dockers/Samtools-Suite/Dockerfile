FROM ubuntu:24.04
RUN apt-get update && apt-get install -y git python3 python3-pandas wget \
libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev curl \
libdeflate-dev bzip2

# Htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
tar -xvf htslib-1.21.tar.bz2 && \
cd htslib-1.21 && \
./configure && \
make && \
make install

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
tar -xvf samtools-1.21.tar.bz2 && \
cd samtools-1.21 && \
./configure && \
make && \
make install

# Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
tar -xvf bcftools-1.21.tar.bz2 && \
cd bcftools-1.21 && \
./configure && \
make && \
make install

RUN curl -sSL https://sdk.cloud.google.com | bash
