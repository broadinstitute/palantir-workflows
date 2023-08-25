FROM ubuntu:22.10

ADD https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary /root/bedtools
RUN chmod 755 /root/bedtools
ENV PATH="${PATH}:/root/"

