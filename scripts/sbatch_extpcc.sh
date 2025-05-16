#!/bin/bash
#SBATCH -D /home/ajmuhich/kliebengrp/ext_pcc/ext_pcc_1/
#SBATCH -o /home/ajmuhich/slurm-log/extpcc_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/extpcc_stderr-%j.txt
#SBATCH -J ext_pcc
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R
R CMD BATCH ~/kliebengrp/ext_pcc/ext_pcc_1/scripts/extract_pcc.R
