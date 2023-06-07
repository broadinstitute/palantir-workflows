FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install apt-utils
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install autoconf
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install libssl-dev

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin
ENV LIB=/usr/local/lib
ENV INC=/usr/local/include

WORKDIR $SRC

# htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
	tar xvf htslib-1.17.tar.bz2 && \
	cd htslib-1.17 && \
	make && \
	make install

# libBigWig
RUN wget https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.7.tar.gz && \
	tar xvf 0.4.7.tar.gz && \
	cd libBigWig-0.4.7 && \
	make && \
	make install

# MethylDackel
RUN wget https://github.com/dpryan79/MethylDackel/archive/refs/tags/0.5.3.tar.gz && \
	tar xvf 0.5.3.tar.gz && \
	cd MethylDackel-0.5.3 && \
	make LIBBIGWIG="/usr/local/lib/libBigWig.a" && \
	make install CFLAGS="-O3 -Wall -I/usr/local/include" LIBS="-L/usr/local/lib" prefix=/usr/local/bin LIBBIGWIG="/usr/local/lib/libBigWig.a"

RUN export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
