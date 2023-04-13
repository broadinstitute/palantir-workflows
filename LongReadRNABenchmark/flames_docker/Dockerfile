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
	python \
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

# Minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2 && \
   tar xvf minimap2-2.24.tar.bz2 && \
   cd minimap2-2.24 && \
   make && \
   cp minimap2 $BIN/

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install

# Pysam
RUN pip install pysam

# numpy
RUN pip install numpy

# editdistance
RUN pip install editdistance

# FLAMES
WORKDIR $SRC

RUN git clone https://github.com/LuyiTian/FLAMES.git
