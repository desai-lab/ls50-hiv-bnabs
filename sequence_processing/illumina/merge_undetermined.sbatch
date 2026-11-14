#!/bin/bash
#SBATCH -J merge_undetermined       # Job name
#SBATCH -c 1                        # Number of cores
#SBATCH -t 0-12:00                   # Runtime (D-HH:MM) - 30 min should be plenty
#SBATCH -p desai           # Partition
#SBATCH --mem=8GB                  
#SBATCH -o ../logs/merge_undetermined_%j.out  # Standard output (%j = job ID)
#SBATCH -e ../logs/merge_undetermined_%j.err  # Standard error

# ---- Configuration ----
# Set these to the directory containing your *_init_counts.csv files
# and the desired output directory for merged files.
INPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/grab_barcodes"
OUTPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_barcodes"

mkdir -p ${OUTPUT_DIR}

# ---- Setup ----
module load python

# ---- Run ----
python merge_undetermined.py "${INPUT_DIR}" "${OUTPUT_DIR}"

