FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install apt-utils
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install unzip
RUN apt-get -qq -y install default-jdk

RUN java -version
RUN javac -version

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
	unzip fastqc_v0.11.9.zip && \
	cd FastQC && \
	chmod +x ./fastqc
