#!/bin/bash
#SBATCH --job-name=extract_bc
#SBATCH --partition=desai
#SBATCH --time=0-00:30:00           # 30 minutes (should be plenty for reading a CSV)
#SBATCH --mem=20G                   # 16GB RAM (adjust if CSV is very large)
#SBATCH --cpus-per-task=1           # Single core is fine for this
#SBATCH --output=../logs/extract_bc_%j.out
#SBATCH --error=../logs/extract_bc_%j.err

#==============================================================================
# EXTRACT BARCODES FROM CSV
# Prepares input file for starcode
#==============================================================================

set -e

echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo ""

# Load Python module (includes mamba)
module load python/3.12.5-fasrc01

# Option 1: If you have an existing conda/mamba environment with pandas
# Uncomment and edit the line below:
# mamba activate /n/holylabs/LABS/<PI_LAB>/Lab/conda/<ENV_NAME>

# Option 2: Use a simple pip install in a temp environment
# Create environment if it doesn't exist
if [ ! -d "$HOME/.conda/envs/barcode_extr_env" ]; then
    echo "Creating conda environment with pandas..."
    mamba create -n barcode_extr_env python=3.11 pandas -y
fi

# Activate the environment
source activate barcode_extr_env

# Verify pandas is available
python -c "import pandas; print(f'pandas version: {pandas.__version__}')"

# Run the extraction script
python extract_barcodes.py

# Deactivate environment
mamba deactivate

echo ""
echo "Finished: $(date)"

