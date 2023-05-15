FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update
RUN apt-get -qq -y install automake
RUN apt-get -qq -y install build-essential
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

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar xvf samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
	./configure && \
	make && \
	make install

# BioPython
RUN pip3 install biopython

# bx-python
RUN pip3 install bx-python

# cython
RUN pip3 install Cython

# numpy should be <1.20, otherwise cDNA_Cupcake fails
RUN pip3 uninstall --yes numpy
RUN pip3 install --upgrade numpy==1.19

# cDNA_Cupcake
WORKDIR $SRC

COPY remove_fastq_duplicates.py $SRC/remove_fastq_duplicates.py

RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git && \
	cd cDNA_Cupcake && \
	python3 setup.py build && \
	python3 setup.py install
