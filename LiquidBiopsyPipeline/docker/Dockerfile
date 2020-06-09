FROM broadinstitute/gatk:4.1.4.1
RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y && \ 
        apt-get install -y python-pandas && \
        apt-get install -y python-biopython && \
        apt-get install -y python-pysam && \
        apt-get clean

#Set the working directory to be used when the docker gets run
WORKDIR /usr
# Install the ggplot2 library and a few other dependencies we want to have available
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('tidyverse')"

# Install HTSLIB
RUN wget 'https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2' -O htslib-1.8.tar.bz2 && tar xf htslib-1.8.tar.bz2 && rm htslib-1.8.tar.bz2 && cd htslib-1.8 && make && make install && cd ../ && rm -r htslib-1.8
RUN wget 'https://github.com/lh3/bwa/releases/download/v0.7.15/bwakit-0.7.15_x64-linux.tar.bz2' -O bwakit-0.7.15.tar.bz2 &&   tar xf bwakit-0.7.15.tar.bz2 &&   rm bwakit-0.7.15.tar.bz2 &&   mv bwa.kit/bwa bwa &&   rm -rf bwa.kit
# Install BLAST 
RUN wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz' -O ncbi-blast-2.2.30+-x64-linux.tar.gz && tar xf ncbi-blast-2.2.30+-x64-linux.tar.gz && rm ncbi-blast-2.2.30+-x64-linux.tar.gz && mv ncbi-blast-2.2.30+/bin/blastn blastn_2.2.30+ && rm -rf ncbi-blast-2.2.30+

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar xf bcftools-1.9.tar.bz2 && rm bcftools-1.9.tar.bz2 && cd bcftools-1.9 && make && make install && cd ../ && rm -r bcftools-1.9

RUN wget 'https://github.com/fulcrumgenomics/fgbio/releases/download/1.0.0/fgbio-1.0.0.jar'

RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" -O GenomeAnalysisTK-3.8-1-0.tar.bz2 && tar xvf GenomeAnalysisTK-3.8-1-0.tar.bz2 && mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar . && rm -r GenomeAnalysisTK-3.8-1-0-gf15c1c3ef

RUN wget 'https://github.com/broadinstitute/picard/releases/download/2.20.8/picard.jar' 
RUN wget 'https://github.com/broadinstitute/picard/releases/download/2.20.8/picardcloud.jar'

COPY . .
