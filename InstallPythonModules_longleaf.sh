#!/bin/bash

module load python/3.5.1
pip3 install --user snakemake multiqc

module load python/2.7.12
pip install --user pysam pyBigWig
