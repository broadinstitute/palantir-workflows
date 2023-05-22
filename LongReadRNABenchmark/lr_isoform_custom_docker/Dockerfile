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

COPY split_gtf.py $SRC/split_gtf.py
COPY summarize_results.py $SRC/summarize_results.py
COPY extract_denovo_model_stats.py $SRC/extract_denovo_model_stats.py
COPY plot_summary_results.py $SRC/plot_summary_results.py
COPY plot_denovo_stats.py $SRC/plot_denovo_stats.py

WORKDIR $SRC
