# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

RUN apt-get update \
  && apt-get install -y wget

RUN wget https://github.com/NovembreLab/admixture/raw/master/releases/admixture_linux-1.3.0.tgz \
  && tar xvzf admixture_linux-1.3.0.tgz \
  && rm admixture_linux-1.3.0.tgz

