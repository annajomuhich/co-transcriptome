#!/bin/bash
#SBATCH -D /home/ajmuhich/bcin_genemodel
#SBATCH -o /home/ajmuhich/slurm-log/bcin_genemodel_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/bcin_genemodel_stderr-%j.txt
#SBATCH -J bcin_genemodel
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

R CMD BATCH scripts/bcin_genemodel.R

