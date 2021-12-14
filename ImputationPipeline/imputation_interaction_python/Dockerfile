# default ubuntu version: 20.04
ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION}

# the following argument is needed for pkg-config to
# install without hanging on user input
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3 \
  && apt-get install -y python3-pip


RUN pip install --upgrade pip
RUN pip install cyvcf2 pandas pandas-plink