#!/bin/bash

source ~/Miniconda3/etc/profile.d/conda.sh

conda activate snakemake

snakemake -n -s Snakefile_gwas_kmer_partI --configfile $1 --use-conda --latency-wait 60 --rerun-incomplete -p -j 15
