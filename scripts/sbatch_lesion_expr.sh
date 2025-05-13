#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/lesion_expr_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/lesion_expr_stderr-%j.txt
#SBATCH -J lesion_expr
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
LESIONSIZE_FILE=$1
BCCOUNTS_FILE=$2
HOST_EXPR_FILE=$3
OUTPUT_DIR=$4

Rscript ~/co-transcriptome/scripts/lesion_expr_model.R \
"$LESIONSIZE_FILE" \
"$BCCOUNTS_FILE" \
"$HOST_EXPR_FILE" \
"$OUTPUT_DIR"