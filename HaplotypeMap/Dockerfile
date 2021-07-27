# 1. Base image
FROM python:3.8.3-slim-buster

# 2. Copy files
COPY . /src

# 3. Install dependencies
RUN pip install argparse
RUN pip install PyVCF
RUN pip install pandas