FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install bc
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install libncurses-dev
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libcurl4-openssl-dev

RUN pip install matplotlib
RUN pip install numpy
RUN pip install pandas
RUN pip install seaborn

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

COPY mapping_efficiency_statistics.sh $SRC/mapping_efficiency_statistics.sh

WORKDIR $SRC

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar xvf samtools-1.10.tar.bz2 && \
	cd samtools-1.10 && \
	make && \
	cp samtools $BIN/
