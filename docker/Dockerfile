FROM ubuntu:24.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-venv \
    python3-dev \
    python3-pip \
    python3-pysam \
    git \
    samtools \
    bcftools \
    libssl-dev \
    vim-common \
    gawk \
    parallel \
    bowtie2 \
    && rm -rf /var/lib/apt/lists/*

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    cd STAR-2.7.11b/source && \
    make STAR && \
    cp STAR /usr/local/bin

# Install seqtk
RUN wget https://github.com/lh3/seqtk/archive/v1.3.tar.gz && \
    tar -xzf v1.3.tar.gz && \
    cd seqtk-1.3 && \
    make && \
    cp seqtk /usr/local/bin

# Set up entry point
CMD ["/bin/bash", "bin/viral_detection.sh", "config/pipeline_input.txt"]
