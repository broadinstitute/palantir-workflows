FROM python:3.11.4

# Install dependencies
RUN pip install numpy==1.25.2
RUN pip install pandas==1.5.3
RUN pip install matplotlib==3.7.1
RUN pip install shapely==1.8.5
RUN pip install descartes==1.1.0

# Install alphashape
RUN pip install alphashape==1.3.1

# Copy custom pca novelty detection script(s)
ENV SRC=/usr/local/src
WORKDIR $SRC
COPY generate_alphashape.py $SRC/generate_alphashape.py
COPY pca_novelty_detection.py $SRC/pca_novelty_detection.py
