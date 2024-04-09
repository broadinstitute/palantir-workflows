# just getting this to work
FROM us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124

RUN apt-get update \
  && apt-get install -y python3 \
  && apt-get install -y python3-pip \
  && apt-get install -y libbz2-dev \
  && apt-get install -y liblzma-dev

RUN pip install pysam==0.20.0
RUN pip install sortedcontainers
RUN pip install pandas

WORKDIR /