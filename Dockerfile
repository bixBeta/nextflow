FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.source="https://github.com/bixBeta/nextflow"
LABEL org.opencontainers.image.description="TREX-RNA: fastp=0.21.0 | STAR=2.7.0e | samtools=1.9 | RSeQC=5.0.1 | fastq_screen=0.15.3 | multiqc=1.32"

USER root

# Layer 1 — large, stable alignment + BAM tools (slowest to download; rarely version-bumped)
RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda \
        star=2.7.0e \
        samtools=1.9 \
    && micromamba clean --all --yes

# Layer 2 — QC / preprocessing tools
RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda \
        fastp=0.21.0 \
        fastq-screen=0.15.3 \
        rseqc=5.0.1 \
    && micromamba clean --all --yes

# Layer 3 — reporting (most likely to be updated independently)
RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda \
        multiqc=1.32 \
    && micromamba clean --all --yes

# --- ADD NEW TOOLS BELOW THIS LINE ---
# Each new RUN block becomes its own cached layer.
# Example:
#   RUN micromamba install -y -n base -c conda-forge -c bioconda \
#           <new-tool>=<version> \
#       && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"
