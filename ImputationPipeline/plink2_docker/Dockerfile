# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

RUN apt-get update \
  && apt-get install -y wget \
  && apt-get install -y unzip

RUN wget https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip
RUN unzip plink2_linux_avx2.zip -d /
RUN rm plink2_linux_avx2.zip

WORKDIR /