FROM ubuntu:24.04
RUN apt-get update && apt-get install -y git python3 python-is-python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev python3-pandas curl

# Make bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
RUN tar -xvf bcftools-1.21.tar.bz2
WORKDIR bcftools-1.21
RUN ./configure && make && make install
RUN export PATH=/bcftools-1.21/bcftools:$PATH

# Install whatshap
RUN pip install whatshap==2.3 --break-system-packages

RUN export PATH=$HOME/.local/bin:$PATH
WORKDIR /
