FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install build-essential
RUN apt-get -qq -y install python3
RUN apt-get -qq -y install python3-dev
RUN apt-get -qq -y install python3-distutils
RUN apt-get -qq -y install python3-pip
RUN apt-get -qq -y install python3-setuptools
RUN apt-get -qq -y install wget

RUN pip install matplotlib
RUN pip install numpy
RUN pip install pandas
RUN pip install seaborn

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

COPY summarize_analysis.py $SRC/summarize_analysis.py
COPY summarize_reffree_analysis.py $SRC/summarize_reffree_analysis.py
COPY summarize_denovo_analysis.py $SRC/summarize_denovo_analysis.py
COPY plot_analysis_summary.py $SRC/plot_analysis_summary.py
COPY plot_denovo_analysis_summary.py $SRC/plot_denovo_analysis_summary.py

WORKDIR $SRC
