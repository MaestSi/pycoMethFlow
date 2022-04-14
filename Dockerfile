FROM continuumio/miniconda3

RUN apt-get update -qq  && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    samtools \
    minimap2 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN conda install python=3.6 pip
RUN python -m pip install pycometh==v0.4.25
RUN conda install -c bioconda nanopolish

WORKDIR /home
