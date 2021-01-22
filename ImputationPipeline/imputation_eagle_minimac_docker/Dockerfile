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


RUN pip3 install docopt
