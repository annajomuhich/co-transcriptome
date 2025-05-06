#!/bin/bash
#SBATCH -D /home/ajmuhich/mr2mods
#SBATCH -o /home/ajmuhich/slurm-log/mr2mods_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/mr2mods_stderr-%j.txt
#SBATCH -J mr2mods
#SBATCH -t 8-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

# ========== USER VARIABLES ==========
INPUT_MATRIX=~/kliebengrp/pvbc_coex/normalized.matrix
OUTPUT_DIR=~/kliebengrp/pvbc_coex
SCRIPT_DIR=/home/ajmuhich/mr2mods/scripts
THREADS=20
D_VALUES=(100 50 25 10 5)
# ====================================

echo "[$(date)] Starting PCC/MR table generation..."
python3 "$SCRIPT_DIR/calc_mr_pcc.py" -i "$INPUT_MATRIX" -o "$OUTPUT_DIR/expr.mr" -t "$THREADS"
echo "[$(date)] Finished PCC/MR table generation."

# Load necessary modules
module load jdk/23.0.1

# Call networks with varying decay rates
cd "$OUTPUT_DIR" || exit 1
for d in "${D_VALUES[@]}"; do
  echo "[$(date)] Starting network/module creation with decay rate $d..."
  python3 "$SCRIPT_DIR/create_network_and_modules.py" -i "$OUTPUT_DIR/expr.mr" -c "$SCRIPT_DIR/cluster_one-1.0.jar" -d "$d"
  echo "[$(date)] Finished network/module creation with decay rate $d."
done

echo "[$(date)] All tasks completed."
