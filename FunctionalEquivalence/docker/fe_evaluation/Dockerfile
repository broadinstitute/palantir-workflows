FROM continuumio/miniconda3

ADD environment.yml /fe_evaluation/environment.yml

RUN conda env create -f /fe_evaluation/environment.yml
RUN echo "source activate fe_evaluation" >> ~/.bashrc