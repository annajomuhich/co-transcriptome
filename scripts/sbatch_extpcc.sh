#!/bin/bash
#SBATCH -D /home/ajmuhich/kliebengrp/ext_pcc/
#SBATCH -o /home/ajmuhich/slurm-log/extpcc_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/extpcc_stderr-%j.txt
#SBATCH -J ext_pcc
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
CLUST_FILE=$1
PCC_DIR=$2
OUTPUT_DIR=$3

Rscript ~/co-transcriptome/scripts/extract_pcc.R \
"$CLUST_FILE" \
"$PCC_DIR" \
"$OUTPUT_DIR"

