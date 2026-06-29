FROM temp_glimpse_base

WORKDIR /docker_build/

# Install required packages
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y openjdk-8-jre-headless curl python3

# Download and build bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
tar -xf bcftools-1.16.tar.bz2 && \
rm bcftools-1.16.tar.bz2 && \
cd bcftools-1.16 && \
autoheader && \
autoconf && \
./configure --enable-libcurl && \
make install && \
cd .. && \
rm -r bcftools-1.16

# Download picard (for updating sequence dictionary)
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && \
mv picard.jar /

WORKDIR /
# Install gcloud
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz && \
tar -xzf google-cloud-cli-linux-x86_64.tar.gz && \
./google-cloud-sdk/install.sh -q --path-update True && \
rm google-cloud-cli-linux-x86_64.tar.gz

RUN rm -r /docker_build/
