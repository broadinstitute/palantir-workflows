FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install apt-utils
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# seqtk
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz && \
	tar xvf v1.3.tar.gz && \
	cd seqtk-1.3 && \
	make && \
	cp seqtk $BIN/
