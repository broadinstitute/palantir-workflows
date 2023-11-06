FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install bzip2
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

# Install alphashape python dependencies
RUN pip install numpy
RUN pip install pandas
RUN pip install matplotlib
RUN pip install shapely==1.8.5
RUN pip install descartes==1.1.0

# Install alphashape
RUN pip install alphashape

# Copy custom pca novelty detection script
COPY pca_novelty_detection.py $SRC/pca_novelty_detection.py
