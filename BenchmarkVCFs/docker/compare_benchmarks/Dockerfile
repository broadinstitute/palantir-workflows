FROM continuumio/miniconda3

ADD environment.yml /compare_benchmarks/environment.yml

RUN conda env create -f /compare_benchmarks/environment.yml
RUN echo "source activate compare_benchmarks" >> ~/.bashrc