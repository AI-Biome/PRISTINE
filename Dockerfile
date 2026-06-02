FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies

RUN apt-get update && apt-get install -y \
    procps \
    curl \
    git \
    build-essential \
    zlib1g-dev \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Conda environment

COPY my_environment.yml /tmp/my_environment.yml

SHELL ["/bin/bash", "-c"]

RUN conda update -n base -c defaults conda -y && \
    conda env create -f /tmp/my_environment.yml && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate primer_design_tool_full && \
    conda install -y pyyaml && \
    conda install -y -c conda-forge -c bioconda pblat && \
    conda clean -afy

# Install OBITools4

RUN cd /opt && \
    curl -L https://raw.githubusercontent.com/metabarcoding/obitools4/master/install_obitools.sh \
    | bash

RUN cd /opt && \
    git clone https://github.com/metabarcoding/obitools4.git && \
    cd obitools4 && \
    sed -i 's/sudo //g' install_obitools.sh && \
    bash install_obitools.sh

# Copy PRISTINE

COPY PRISTINE_essentials.py /opt/PRISTINE_essentials.py

RUN chmod +x /opt/PRISTINE_essentials.py

# Environment

ENV PATH=/opt/conda/envs/primer_design_tool_full/bin:$PATH
ENV CONDA_DEFAULT_ENV=primer_design_tool_full

# Default execution

ENTRYPOINT ["/bin/bash", "-c", "\
source /opt/conda/etc/profile.d/conda.sh && \
conda activate primer_design_tool_full && \
exec python /opt/PRISTINE_essentials.py \"$@\"", "--"]
