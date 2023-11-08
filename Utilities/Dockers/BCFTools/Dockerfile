FROM ubuntu:23.10
RUN apt-get update && apt-get install -y git python3 python-is-python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev python3-pandas

# Make bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
RUN tar -xvf bcftools-1.18.tar.bz2
WORKDIR bcftools-1.18
RUN ./configure && make && make install
RUN export PATH=/bcftools-1.18/bcftools:$PATH
