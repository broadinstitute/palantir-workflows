FROM ubuntu:22.10
RUN apt update
RUN apt install -y git python3 python3-pip wget libbz2-dev liblzma-dev libcurl4-openssl-dev

RUN pip install pandas

# Make bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
RUN tar -xvf bcftools-1.16.tar.bz2
WORKDIR bcftools-1.16
RUN ./configure
RUN make
RUN make install
RUN export PATH=/bcftools-1.16/bcftools:$PATH

RUN curl -sSL https://sdk.cloud.google.com | bash