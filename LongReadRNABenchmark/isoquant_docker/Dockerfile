FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libncurses-dev
RUN apt-get -qq -y install libncurses5-dev
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2 && \
	tar xvf minimap2-2.24.tar.bz2 && \
	cd minimap2-2.24 && \
	make && \
	cp minimap2 $BIN/

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
	./configure && \
	make && \
	make install

# IsoQuant 3.2.0 for memory optimizations
RUN wget https://github.com/ablab/IsoQuant/releases/download/v3.2.0/IsoQuant-3.2.0.tar.gz && \
	tar xvf IsoQuant-3.2.0.tar.gz && \
	cd IsoQuant-3.2.0 && \
	pip install numpy==1.21 && \
	pip install -r requirements.txt
