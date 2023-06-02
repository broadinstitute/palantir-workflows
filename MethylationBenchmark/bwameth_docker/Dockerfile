FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install apt-utils
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install libncurses-dev
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install unzip

# toolshed
RUN pip3 install toolshed==0.4.6

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# bwa
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
   tar xvf bwa-0.7.17.tar.bz2 && \
   cd bwa-0.7.17 && \
   make && \
   cp bwa $BIN/

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar xvf samtools-1.10.tar.bz2 && \
	cd samtools-1.10 && \
	make && \
	cp samtools $BIN/

# bwameth
RUN wget https://github.com/brentp/bwa-meth/archive/refs/tags/v0.2.2.zip && \
	unzip v0.2.2.zip && \
	cd bwa-meth-0.2.2 && \
	python3 setup.py install
