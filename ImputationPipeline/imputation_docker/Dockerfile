FROM ubuntu:latest

RUN apt-get update \
  && apt-get install -y wget \
  && apt-get install -y libgomp1 \
  && apt-get install -y build-essential \
  && apt-get install -y libz-dev \
  && apt-get install -y libbz2-dev \
  && apt-get install -y liblzma-dev \
  && apt-get install -y pkg-config \
  && apt-get install -y tabix \
  && apt-get install -y python3 \
  && apt-get install -y python3-pip \
  && rm -rf /var/lib/apt/lists/*


RUN wget https://github.com/genepi/imputationserver/raw/master/files/bin/eagle
RUN wget https://github.com/genepi/imputationserver/raw/master/files/bin/Minimac4
RUN ["chmod", "+x", "eagle"]
RUN ["chmod", "+x", "Minimac4"]

RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
  && tar xf bcftools-1.10.2.tar.bz2 \
  && cd bcftools-1.10.2 \
  && ./configure \
  && make \ 
  && make install \
  && cd ../ \
  && rm -r bcftools-1.10.2.tar.bz2

RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz

RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
  && tar xf vcftools-0.1.16.tar.gz \
  && cd vcftools-0.1.16 \
  && ./configure \
  && make \
  && make install \
  && cd ../ \
  && rm -r vcftools-0.1.16.tar.gz

RUN pip3 install docopt
