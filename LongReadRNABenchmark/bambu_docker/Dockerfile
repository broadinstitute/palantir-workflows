FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential 
RUN apt-get -qq -y install curl
RUN apt-get -qq -y install gdebi
RUN apt-get -qq -y install libcurl4-openssl-dev
RUN apt-get -qq -y install wget

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
