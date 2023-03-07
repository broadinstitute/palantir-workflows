FROM ubuntu:22.10

RUN apt update
RUN apt install -y wget
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary bedtools
RUN chmod a+x bedtools
RUN mv bedtools /root/bedtools
ENV PATH="${PATH}:/root/"

