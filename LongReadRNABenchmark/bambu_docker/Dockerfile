FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    ftp \
    g++ \
    gcc \
	gdebi \
	gdebi-core \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev


ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

ARG R_VERSION="4.2.3"

RUN curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb

# --non-interactive == DANGER_ZONE
# Should be fine here though
RUN gdebi --non-interactive r-${R_VERSION}_1_amd64.deb

RUN cp /opt/R/${R_VERSION}/bin/* $BIN/

RUN R -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'

RUN R -e 'BiocManager::install("bambu")'
