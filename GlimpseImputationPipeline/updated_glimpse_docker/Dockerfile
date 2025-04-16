FROM us.gcr.io/broad-dsde-methods/glimpse:odelaneau_bd93ade

# Update bcftools to latest
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
tar -xf bcftools-1.21.tar.bz2 && \
rm bcftools-1.21.tar.bz2 && \
cd bcftools-1.21 && \
autoheader && \
autoconf && \
./configure --enable-libcurl && \
make install && \
cd .. && \
rm -r bcftools-1.21

# Install gcloud
# RUN curl -sSL https://sdk.cloud.google.com | bash

# Make gcloud commands available
ENV PATH=$PATH:/root/google-cloud-sdk/bin