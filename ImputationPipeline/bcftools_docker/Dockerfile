ARG ALPINE_VERSION=3.12

FROM alpine:${ALPINE_VERSION} AS alpine-build-bcftools
RUN echo http://dl-cdn.alpinelinux.org/alpine/edge/main >> /etc/apk/repositories
RUN echo http://dl-cdn.alpinelinux.org/alpine/edge/community >> /etc/apk/repositories

ARG ALPINE_VERSION=3.12
ARG BCFTOOLS_VERSION=1.10.2

RUN apk update && apk add build-base wget zlib-dev tar bzip2-dev xz-dev \
    curl-dev curl-static zlib-static bzip2-static openssl-libs-static  \
    nghttp2-static brotli-static

WORKDIR /build
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2
WORKDIR /build/bcftools-${BCFTOOLS_VERSION}
RUN ./configure && \
    make CFLAGS="-g -O2 -static" \
    LIBS="-llzma -lssl -lcrypto -lnghttp2 -lcurl -lbrotlienc-static -lbrotlidec-static -lbrotlicommon-static" \
    LDFLAGS="-static"
RUN strip bcftools

FROM alpine:${ALPINE_VERSION}

ARG BCFTOOLS_VERSION=1.10.2

WORKDIR / 
COPY --from=alpine-build-bcftools /build/bcftools-${BCFTOOLS_VERSION}/bcftools /usr/local/bin/bcftools
ENTRYPOINT ["/usr/local/bin/bcftools"]