#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/bcin_genemodel_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/bcin_genemodel_stderr-%j.txt
#SBATCH -J bcin_genemodel
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
COUNTS_FILE1=$1
COUNTS_FILE2=$2
SAMPLEID_FILE1=$3
SAMPLEID_FILE2=$4
BATCH_FILE1=$5
BATCH_FILE2=$6
OUTPUT_DIR=$7

Rscript ~/co-transcriptome/scripts/bcin_genemodel.R \
"$COUNTS_FILE1" "$COUNTS_FILE2" \
"$SAMPLEID_FILE1" "$SAMPLEID_FILE2" \
"$BATCH_FILE1" "$BATCH_FILE2" \
"$OUTPUT_DIR"

