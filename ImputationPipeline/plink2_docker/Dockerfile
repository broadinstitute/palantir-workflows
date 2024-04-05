# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

RUN apt-get update \
  && apt-get install -y wget \
  && apt-get install -y unzip \
  && apt-get install -y python3 \
  && apt-get install -y python3-pip \
  && apt-get install -y libbz2-dev \
  && apt-get install -y liblzma-dev

RUN pip install pysam==0.20.0
RUN pip install sortedcontainers
RUN pip install pandas

RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip
RUN unzip plink2_linux_avx2_20240105.zip -d /
RUN rm plink2_linux_avx2_20240105.zip

WORKDIR /