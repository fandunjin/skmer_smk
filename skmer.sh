#!/bin/bash
#JSUB -J skmer
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_log.%J
#JSUB -e skmer_err.%J
#JSUB -cwd skmer_%J

source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
conda activate 01bio

conda env list > condainfo

cd /hpcfile/users/92024286/skmer_smk/

snakemake -s /hpcfile/users/92024286/skmer_smk/snakefile --core 48 --rerun-incomplete --latency-wait 120
