FROM ubuntu:22.04
LABEL author=dilldesilva02@gmail.com
#install all necessary things here
RUN apt-get -y update
RUN apt-get -y install curl wget gcc libz-dev liblzma-dev libbz2-dev libcurl4-gnutls-dev libhts-dev libhts-dev git g++ gzip bzip2 make
#set-up of htslib-1.17
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
RUN tar -xvf htslib-1.17.tar.bz2
WORKDIR /htslib-1.17
RUN ./configure --prefix=/usr/local 
RUN make
RUN make install
RUN export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
#set-up of vcfdist
WORKDIR /
RUN git clone https://github.com/timd1/vcfdist
#RUN cd htslib-1.17/vcfdist/src
WORKDIR /vcfdist/src
RUN make
RUN cp /vcfdist/src/vcfdist /usr/local/bin
WORKDIR /