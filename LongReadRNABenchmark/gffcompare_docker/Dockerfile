FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install zlib1g-dev

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# GFFCompare
RUN wget https://github.com/gpertea/gffcompare/releases/download/v0.12.6/gffcompare-0.12.6.tar.gz && \
	tar xvf gffcompare-0.12.6.tar.gz && \
	cd gffcompare-0.12.6 && \
	make release && \
	cp gffcompare $BIN/
