#!/bin/bash

source ~/Miniconda3/etc/profile.d/conda.sh

conda activate snakemake

snakemake -s Snakefile_gwas_kmer_partII --configfile $1 --use-conda --latency-wait 60 --rerun-incomplete -p -j 15

conda deactivate

echo 'PART II FINISHED SUCCESSFULLY'
