# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

ARG PLINK_VERSION=v2.00a2.3

RUN apt-get update \
  && apt-get install -y wget unzip build-essential zlib1g-dev libatlas-base-dev

RUN wget https://github.com/chrchang/plink-ng/archive/refs/tags/$PLINK_VERSION.zip
RUN unzip $PLINK_VERSION.zip
RUN plink_dir=$(find . -type d -name "plink*") && \
    make -C ${plink_dir}/2.0/build_dynamic/ NO_AVX2= && \
    mv ${plink_dir}/2.0/build_dynamic/plink2 / && \
    rm -r ${plink_dir}

WORKDIR /