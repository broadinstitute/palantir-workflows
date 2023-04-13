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

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install

# BioPython
RUN pip3 install biopython

# bx-python
RUN pip3 install bx-python

# cython
RUN pip3 install Cython

# setuptools
RUN apt-get install python3-setuptools

# numpy should be <1.20, otherwise cDNA_Cupcake fails
RUN pip3 uninstall --yes numpy
RUN pip3 install --upgrade numpy==1.19

# cDNA_Cupcake
WORKDIR $SRC

COPY remove_fastq_duplicates.py $SRC/remove_fastq_duplicates.py

RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git && \
	cd cDNA_Cupcake && \
	python3 setup.py build && \
	python3 setup.py install
