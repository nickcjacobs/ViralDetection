# 🦠 Viral Detection Pipeline

This repository contains a reproducible pipeline for detecting viral sequences from sequencing data.  
It can be run using **Docker** (recommended for full reproducibility) or alternatively with **Conda** if Docker is unavailable.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Option 1: Run with Docker (Recommended)](#-option-1-run-with-docker-recommended)
- [Option 2: Run with Conda](#-option-2-run-with-conda)
- [Input Configuration](#-input-configuration)
- [Output](#-output)
- [Updating or Removing Environments](#-updating-or-removing-environments)

---

## 🧬 Overview

This pipeline performs the following key steps:
- Quality control and alignment of sequencing reads
- Viral sequence identification using `bowtie2` and `STAR`
- Variant calling with `bcftools`
- Optional downstream analysis and summary statistics

The codebase supports:
- **Docker** for fully containerized execution  
- **Conda** for systems where Docker cannot be used

---

## 🐋 Option 1: Run with Docker (Recommended)

Docker ensures complete reproducibility with all dependencies pre-installed.

### Prerequisites
- [Docker installed](https://docs.docker.com/get-docker/)
- Internet connection for image building (first time only)

### Clone this repository

```
git clone https://github.com/nickcjacobs/ViralDetection.git
cd ViralDetection
```

### Build the Docker image

From the repository root:

```
docker build -t viral_detection .
```

### Run the pipeline

```
docker run --rm -v $(pwd):/app -w /app viral_detection \
    bash bin/viral_detection.sh config/pipeline_input.txt
```

Explanation:

-v $(pwd):/app mounts your current directory into the container

-w /app sets the working directory inside the container

The pipeline reads parameters from config/pipeline_input.txt

## 🧫 Option 2: Run with Conda

If Docker is not available, you can run the same pipeline in a Conda environment.

### Prerequisites

Miniconda
 or
Mambaforge

### Step 1: Clone this repository
```
git clone https://github.com/nickcjacobs/ViralDetection.git
cd ViralDetection
```

### Step 2: Create and activate the environment
```
conda env create -f environment.yml
conda activate viral_detection
```

This installs:

samtools, bcftools, bowtie2, STAR, seqtk, parallel, pysam, and other required tools.

### Step 3: Run the pipeline
```
bash bin/viral_detection.sh config/pipeline_input.txt
```

## ⚙️ Input Configuration

All input file paths and settings are specified in:

config/pipeline_input.txt


Ensure this file includes the correct paths to your FASTQ files, reference genome, and other required inputs before running the pipeline.

## 📊 Output

The pipeline produces:

Processed and aligned reads

Detected viral sequences

Variant calls (.vcf files)

Summary and log files in the designated output folder

Output locations and naming conventions are controlled by your configuration file.

## 🔄 Updating or Removing Environments

If you modify environment.yml and want to apply updates:

```
conda env update -f environment.yml --prune
```

To remove the Conda environment entirely:

```
conda remove --name viral_detection --all
```

Pull requests are welcome!
If you’d like to add new features or improve existing ones:

Fork this repository

Create a feature branch

Submit a pull request describing your changes

For major updates, please open an issue first to discuss proposed modifications.

Maintainer: Nick Jacobs
Repository: github.com/KlugerLab/ViralDetection
