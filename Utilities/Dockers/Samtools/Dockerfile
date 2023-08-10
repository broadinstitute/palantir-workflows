FROM ubuntu:22.10
RUN apt-get update && apt-get install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev curl

RUN pip install pandas

# Make Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
RUN tar -xvf samtools-1.16.1.tar.bz2
WORKDIR samtools-1.16.1
RUN ./configure
RUN make
RUN make install
RUN export PATH=/samtools-1.16.1/samtools:$PATH

RUN curl -sSL https://sdk.cloud.google.com | bash
