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

Rscript scripts/clustertopGO.R \
  "$1" \  # host annotation file path
  "$2" \  # bcin annotation file path
  "$3" \  # cluster data file path
  "$4"    # output directory path

