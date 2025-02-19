FROM ubuntu:24.04
RUN apt update
RUN apt install -y openjdk-11-jre ant git unzip python3 python3-pip libbz2-dev liblzma-dev libcurl4-openssl-dev python3-pandas wget

# Get RTG files for install
RUN git clone https://github.com/RealTimeGenomics/rtg-tools.git
WORKDIR rtg-tools
RUN ant zip-nojre
WORKDIR dist
RUN unzip *.zip

# Install RTG
WORKDIR rtg-tools-*
RUN ln -s $(pwd)/rtg /bin/
RUN yes no | rtg --help
WORKDIR /root

# Add bedtools
WORKDIR /
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static .
RUN mv bedtools.static bedtools && chmod +x bedtools && cp bedtools /bin/

# Htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
tar -xvf htslib-1.21.tar.bz2 && \
cd htslib-1.21 && \
./configure && \
make && \
make install

# Make bcftools
ADD https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 .
RUN tar -xvf bcftools-1.21.tar.bz2
WORKDIR bcftools-1.21
RUN ./configure
RUN make
RUN make install
WORKDIR /
