FROM rocker/tidyverse:4.4
ARG BIOCONDUCTOR_VERSION=3.20
RUN R -e 'install.packages("BiocManager")'
RUN R -e "BiocManager::install('rhdf5', version='${BIOCONDUCTOR_VERSION}')"