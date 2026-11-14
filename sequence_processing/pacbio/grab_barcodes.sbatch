#!/bin/bash
#SBATCH -J barcode_extract          # Job name
#SBATCH -p desai           # Partition (serial_requeue is efficient for single-core jobs)
#SBATCH -n 1                        # Number of cores
#SBATCH -t 0-10:00                  # Runtime (D-HH:MM) - adjust based on your file size
#SBATCH --mem=16G                    # Memory allocation - adjust if needed for large FASTQ files
#SBATCH -o ../logs/grab_barcodes_%j.out   # Standard output (%j = job ID)
#SBATCH -e ../logs/grab_barcodes_%j.err   # Standard error (%j = job ID)


# Print job info for debugging
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Job ID: $SLURM_JOB_ID"

# Load the python module
module load python

# Activate the conda environment with your packages
source activate barcode_env

# Verify packages are available (optional - useful for debugging)
echo "Python version: $(python --version)"
echo "Checking required packages..."
python -c "import numpy; import pandas; from Bio.SeqIO.QualityIO import FastqGeneralIterator; import edlib; print('All packages loaded successfully')"

# Run the script
python grab_barcodes.py

# Print completion info
echo "Job completed at $(date)"

