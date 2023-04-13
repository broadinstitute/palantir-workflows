FROM ubuntu:20.04

# TALON
# https://github.com/mortazavilab/TALON

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

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
        tar xvf samtools-1.16.1.tar.bz2 && \
        cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install


RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN pip install cython pytest

RUN  wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
        chmod 755 bedtools.static.binary && \
        mv bedtools.static.binary /usr/local/bin/bedtools

# Talon
RUN wget https://github.com/mortazavilab/TALON/archive/refs/tags/v5.0.tar.gz && \
    tar xvf v5.0.tar.gz && \
    cd TALON-5.0 && \
    pip install .
