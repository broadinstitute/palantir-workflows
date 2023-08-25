FROM ubuntu:22.10
RUN apt update
RUN apt install -y openjdk-11-jre ant git unzip python3 python3-pip

# Get RTG files for install
RUN git clone https://github.com/RealTimeGenomics/rtg-tools.git
WORKDIR rtg-tools
RUN ant zip-nojre
WORKDIR dist
RUN unzip *.zip 

# Install RTG
WORKDIR rtg-tools-*
RUN ln -s $(pwd)/rtg /bin/
RUN yes no | rtg --help
WORKDIR /root

# Install Python packages
RUN pip install pandas
