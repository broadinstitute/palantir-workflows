FROM ubuntu:24.04
RUN apt-get update && apt-get install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev curl

RUN pip install pandas

# Make Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
tar -xvf samtools-1.16.1.tar.bz2 && \
cd samtools-1.16.1 && \
./configure && \
make && \
make install

RUN curl -sSL https://sdk.cloud.google.com | bash
