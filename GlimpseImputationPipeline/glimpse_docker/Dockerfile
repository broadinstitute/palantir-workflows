FROM ubuntu:20.04

LABEL org.opencontainers.image.created="2022-11-30"
LABEL org.opencontainers.image.url="https://github.com/odelaneau/GLIMPSE"
LABEL org.opencontainers.image.version="2.0.0"
LABEL org.opencontainers.image.licences="MIT"
LABEL org.opencontainers.image.title="glimpse"
LABEL org.opencontainers.image.authors="simone.rubinacci@unil.ch"

WORKDIR /docker_build/

# Install required packages
RUN apt-get update && apt-get install -y build-essential libbz2-dev libcurl4-openssl-dev autoconf libssl-dev wget zlib1g-dev liblzma-dev git

# Download and build boost program_options and iostreams
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz && \
tar -xf boost_1_78_0.tar.gz && \
rm boost_1_78_0.tar.gz && \
cd boost_1_78_0/ && \
./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
./b2 install && \
cd .. && \
cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
cp -r boost/include/boost/ /usr/include/ && \
rm -r boost_1_78_0 boost

# Download and build htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
tar -xf htslib-1.16.tar.bz2 && \
rm htslib-1.16.tar.bz2 && \
cd htslib-1.16 && \
autoheader && \
autoconf && \
./configure --enable-libcurl && \
make install && \
cd .. && \
rm -r htslib-1.16

# Download and build GLIMPSE
RUN git clone --depth 1 --branch v2.0.0_fixed_makefiles_for_docker https://github.com/michaelgatzen/GLIMPSE.git && \
cd GLIMPSE && \
make && \
cd .. && \
mkdir /GLIMPSE && \
mv GLIMPSE/chunk/bin/GLIMPSE2_chunk GLIMPSE/split_reference/bin/GLIMPSE2_split_reference GLIMPSE/phase/bin/GLIMPSE2_phase GLIMPSE/ligate/bin/GLIMPSE2_ligate GLIMPSE/concordance/bin/GLIMPSE2_concordance /GLIMPSE && \
rm -r GLIMPSE

# Install java and download picard (for updating sequence dictionary)
RUN apt-get install -y openjdk-8-jre-headless
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && \
mv picard.jar /

WORKDIR /