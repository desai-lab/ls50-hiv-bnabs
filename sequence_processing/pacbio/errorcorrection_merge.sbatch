#!/bin/bash
#SBATCH --job-name=barcode_ec_merge   # Job name
#SBATCH --partition=desai             # Partition
#SBATCH --time=0-01:00                # Time limit (merging is fast)
#SBATCH --mem=16G                     # Memory for loading all chunks
#SBATCH --cpus-per-task=1             # Single threaded
#SBATCH --output=../logs/barcode_ec_merge_%j.out
#SBATCH --error=../logs/barcode_ec_merge_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaoascensao@g.harvard.edu

# ============================================================================
# PARALLEL BARCODE ERROR CORRECTION - MERGE JOB
# ============================================================================
# This script merges results from all array tasks.
# Run this AFTER the array job completes.
#
# Usage:
#   Option 1 (manual): sbatch errorcorrection_merge.sbatch
#   Option 2 (automatic): Use --dependency flag when submitting:
#     ARRAY_JOB_ID=$(sbatch --parsable errorcorrection_array.sbatch)
#     sbatch --dependency=afterok:$ARRAY_JOB_ID errorcorrection_merge.sbatch
# ============================================================================

# Number of tasks (must match array job)
NUM_TASKS=50

# Print job information
echo "========================================"
echo "Merge Job ID: $SLURM_JOB_ID"
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

# Verify environment
echo "Python: $(which python)"
echo ""

# Check that chunk files exist
CHUNK_COUNT=$(ls ../chunks/chunk_*.csv 2>/dev/null | wc -l)
echo "Found $CHUNK_COUNT chunk files in ../chunks/"

if [ "$CHUNK_COUNT" -eq 0 ]; then
    echo "ERROR: No chunk files found! Did the array job complete?"
    exit 1
fi

echo "Starting merge..."
echo "========================================"

# Run the merge script
python errorcorrection_merge.py --chunk-dir ../chunks --num-tasks $NUM_TASKS

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Merge completed successfully!"
    echo "Output files:"
    ls -la ../superbarcodes_all.csv ../superbarcodes_pure.csv ../consensus_seqs.csv 2>/dev/null
    echo ""
    echo "You can now optionally clean up chunk files:"
    echo "  rm -rf ../chunks/"
else
    echo ""
    echo "========================================"
    echo "ERROR: Merge failed with exit code $?"
fi

echo "End time: $(date)"
echo "========================================"
