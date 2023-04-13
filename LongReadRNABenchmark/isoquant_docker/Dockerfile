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

## Minimap2

RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2 && \
   tar xvf minimap2-2.24.tar.bz2 && \
   cd minimap2-2.24 && \
   make && \
   cp minimap2 $BIN/

## Samtools

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install

## IsoQuant

RUN wget https://github.com/ablab/IsoQuant/releases/download/v3.1.1/IsoQuant-3.1.1.tar.gz && \
   tar xvf IsoQuant-3.1.1.tar.gz && \
   cd IsoQuant-3.1.1 && \
   pip install numpy==1.19 && \
   pip install -r requirements.txt






