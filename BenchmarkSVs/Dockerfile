# This Dockerfile is used to build the monolithic docker image used for most of the tasks in the BenchmarkSVs workflow.
# It contains bcftools, truvari, pandas, and bedtools.
FROM us.gcr.io/broad-dsde-methods/bcftools:v1.0
RUN pip install truvari==4.0.0 pandas==1.5.3 pysam==0.21.0

# Add bedtools
WORKDIR /
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static .
RUN mv bedtools.static bedtools && chmod +x bedtools && cp bedtools /bin/

# Add tabix/bgzip
RUN apt-get update && apt-get install -y tabix