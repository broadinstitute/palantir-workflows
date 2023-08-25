FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install bzip2
RUN apt-get -qq -y install cmake
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install ftp
RUN apt-get -qq -y install g++
RUN apt-get -qq -y install gcc
RUN apt-get -qq -y install gdebi
RUN apt-get -qq -y install gdebi-core
RUN apt-get -qq -y install gfortran
RUN apt-get -qq -y install git
RUN apt-get -qq -y install libblas-dev
RUN apt-get -qq -y install libbz2-dev
RUN apt-get -qq -y install libcairo2-dev
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install libdb-dev
RUN apt-get -qq -y install libghc-zlib-dev
RUN apt-get -qq -y install libjpeg-dev
RUN apt-get -qq -y install liblzma-dev
RUN apt-get -qq -y install libncurses-dev
RUN apt-get -qq -y install libncurses5-dev
RUN apt-get -qq -y install libpcre3-dev
RUN apt-get -qq -y install libpng-dev
RUN apt-get -qq -y install libreadline-dev
RUN apt-get -qq -y install libreadline-dev
RUN apt-get -qq -y install libssl-dev
RUN apt-get -qq -y install libtbb-dev
RUN apt-get -qq -y install libx11-dev
RUN apt-get -qq -y install libxml2-dev
RUN apt-get -qq -y install libxt-dev
RUN apt-get -qq -y install libzmq3-dev
RUN apt-get -qq -y install make
RUN apt-get -qq -y install perl
RUN apt-get -qq -y install pkg-config
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install rsync
RUN apt-get -qq -y install unzip
RUN apt-get -qq -y install wget
RUN apt-get -qq -y install x11-common
RUN apt-get -qq -y install zlib1g-dev

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

ARG R_VERSION="4.2.3"

RUN curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb

# --non-interactive == DANGER_ZONE
# Should be fine here though
RUN gdebi --non-interactive r-${R_VERSION}_1_amd64.deb

RUN cp /opt/R/${R_VERSION}/bin/* $BIN/

RUN R -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'

RUN R -e 'BiocManager::install("bambu")'
