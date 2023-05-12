FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install apt-utils
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install zlib1g-dev
RUN apt-get -qq -y install pigz
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install unzip

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# cutadapt
RUN pip3 install cutadapt==3.2

# FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
	unzip fastqc_v0.11.9.zip && \
	cd FastQC && \
	chmod +x ./fastqc && \
	cp fastqc $BIN/

# TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz && \
	tar xvf trim_galore.tar.gz && \
	cp TrimGalore-0.6.5/trim_galore $BIN/
