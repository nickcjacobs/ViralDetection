# My Pipeline

This pipeline processes sequencing data to identify viral integrations.

## Quick Start

```bash
git clone https://github.com/KlugerLab/ViralDetection.git
cd my-pipeline
docker build -f docker/Dockerfile -t viral_detection .
docker run --rm viral_detection
