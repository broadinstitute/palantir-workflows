FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install bzip2
RUN apt-get -qq -y install cmake
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install git
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

# StringTie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz && \
	tar xvf stringtie-2.2.1.tar.gz && \
	cd stringtie-2.2.1/ && \
	make release && \
	cp stringtie $BIN/
