#!/bin/bash
#SBATCH -D /home/ajmuhich/model_means/IT/
#SBATCH -o /home/ajmuhich/slurm-log/model_means_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/model_means_stderr-%j.txt
#SBATCH -J model_means
#SBATCH -t 1-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

#set -x  # Enable debugging (print each command before execution)

module load R

# Define input file paths
COUNTS_FILE=$1
SAMPLEID_FILE=$2
BATCH_FILE=$3
OUTPUT_DIR=$4

Rscript ~/co-transcriptome/scripts/model_means.R "$COUNTS_FILE" "$SAMPLEID_FILE" "$BATCH_FILE" "$OUTPUT_DIR"

