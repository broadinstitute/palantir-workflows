FROM ubuntu:22.10
RUN apt-get update && apt-get install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev curl

# Htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
tar -xvf htslib-1.18.tar.bz2 && \
cd htslib-1.18 && \
./configure && \
make && \
make install

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
tar -xvf samtools-1.18.tar.bz2 && \
cd samtools-1.18 && \
./configure && \
make && \
make install

# Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && \
tar -xvf bcftools-1.18.tar.bz2 && \
cd bcftools-1.18 && \
./configure && \
make && \
make install

RUN curl -sSL https://sdk.cloud.google.com | bash

RUN pip install pandas
