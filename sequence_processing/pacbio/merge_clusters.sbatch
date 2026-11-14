#!/bin/bash
#SBATCH --job-name=barcode_ec         # Job name
#SBATCH --partition=desai            # Partition (shared is typical for single-node jobs)
#SBATCH --time=0-24:00                # Time limit D-HH:MM (adjust based on data size)
#SBATCH --mem=20G                     # Memory allocation (adjust based on data size)
#SBATCH --cpus-per-task=1             # Number of CPU cores (script is single-threaded)
#SBATCH --output=../logs/merge_cl_%j.out    # Standard output file (%j = job ID)
#SBATCH --error=../logs/merge_cl_%j.err     # Standard error file


# Print job information
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $HOSTNAME"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "========================================"
echo ""

# Load required modules
module purge
module load python/3.10.13-fasrc01
module load Mambaforge/23.11.0-fasrc01

# Activate conda environment
source activate ec_env

python merge_clusters.py