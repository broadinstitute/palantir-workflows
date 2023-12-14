FROM ubuntu:23.04
RUN apt-get update && \
    apt-get install -y git python3 python3-pip python3-dev build-essential libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev tabix && \
    pip install whatshap==1.7 && \
    rm -rf /var/lib/apt/lists/*

RUN export PATH=$HOME/.local/bin:$PATH
