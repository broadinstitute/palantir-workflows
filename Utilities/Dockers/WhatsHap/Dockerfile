FROM ubuntu:23.04
RUN apt update
RUN apt install -y git python3 python3-pip python3-dev build-essential libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev tabix

RUN pip install whatshap==1.7

RUN export PATH=$HOME/.local/bin:$PATH
