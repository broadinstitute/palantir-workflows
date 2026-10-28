FROM condaforge/miniforge3:26.9.4-0

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# Install everything in one transaction so mamba never runs after its own dependencies are upgraded.
# The isoseq package also provides the isoseq3 symlink used by IsoSeq.wdl.
RUN mamba install -y -c conda-forge -c bioconda \
        samtools=1.25 \
        pbmm2=26.2.99 \
        isoseq=26.2.0 \
    && mamba clean -afy

RUN samtools --version && pbmm2 --version && isoseq3 collapse --help
