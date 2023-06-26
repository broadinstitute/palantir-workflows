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

RUN pip3 install deeptools

RUN multiBigwigSummary --version
