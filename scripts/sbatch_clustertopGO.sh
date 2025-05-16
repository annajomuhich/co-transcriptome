#!/bin/bash
#SBATCH -D /home/ajmuhich/clustertopGO
#SBATCH -o /home/ajmuhich/slurm-log/clustertopGO_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/clustertopGO_stderr-%j.txt
#SBATCH -J clustertopGO
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
HOSTANNOT_FILE=$1
BCINANNOT_FILE=$2
CLUST_FILE=$3
OUTPUT_DIR=$4

Rscript ~/co-transcriptome/scripts/clustertopGO.R \
"$HOSTANNOT_FILE" \
"$BCINANNOT_FILE" \
"$CLUST_FILE" \
"$OUTPUT_DIR"
