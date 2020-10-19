# Pull base image
FROM alpine:3.12

# Installation instructions: http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200921.zip


RUN apk update && apk add unzip wget build-base && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200921.zip && \
    unzip plink_linux_x86_64_20200921.zip -d /usr/local/ && \
    rm plink_linux_x86_64_20200921.zip && \
    cd /usr/local

# WORKAROUND: plink hangs if not started with '--noweb' option
# See issue #1
RUN echo '#!/bin/sh'                                                 > /usr/local/bin/plink && \
    echo '#Launch the real plink script forcing the --noweb argument' >> /usr/local/bin/plink && \
    echo '/usr/local/plink --noweb "$@"'                        >> /usr/local/bin/plink && \
    chmod a+x /usr/local/bin/plink


ARG BCFTOOLS_VERSION=1.10.2

WORKDIR /build

RUN apk add zlib-dev tar bzip2-dev xz-dev \
    curl-dev curl-static zlib-static bzip2-static openssl-libs-static  \
    nghttp2-static brotli-static
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2
WORKDIR /build/bcftools-${BCFTOOLS_VERSION}
RUN ./configure && \
    make CFLAGS="-g -O2 -static" \
    LIBS="-llzma -lssl -lcrypto -lnghttp2 -lcurl -lbrotlienc-static -lbrotlidec-static -lbrotlicommon-static" \
    LDFLAGS="-static"
RUN strip bcftools
RUN cp bcftools /usr/local/bin && cd /build && rm -r bcftools-${BCFTOOLS_VERSION}
WORKDIR /usr/local/bin

# Set the default action to print plink's options
#ENTRYPOINT ["/usr/local/bin/plink"]