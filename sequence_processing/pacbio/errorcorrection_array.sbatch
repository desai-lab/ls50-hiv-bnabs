#!/bin/bash
#SBATCH --job-name=barcode_ec_array   # Job name
#SBATCH --partition=desai             # Partition
#SBATCH --time=0-12:00                # Time limit per task (reduced since work is split)
#SBATCH --mem=12G                      # Memory per task (reduced since data is shared)
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --array=1-50                  # Array job: 50 parallel tasks (adjust as needed)
#SBATCH --output=../logs/barcode_ec_%A_%a.out   # %A=array job ID, %a=task ID
#SBATCH --error=../logs/barcode_ec_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaoascensao@g.harvard.edu

# ============================================================================
# PARALLEL BARCODE ERROR CORRECTION - ARRAY JOB
# ============================================================================
# This script runs as a SLURM array job. Each task processes a subset of
# super_bcs. After all tasks complete, run the merge script.
#
# To adjust parallelism, modify the --array parameter above.
# Example: --array=1-100 for 100 parallel tasks
#
# Usage:
#   1. Submit this job: sbatch errorcorrection_array.sbatch
#   2. After completion, submit merge job: sbatch errorcorrection_merge.sbatch
# ============================================================================

# Number of tasks (must match --array upper bound)
NUM_TASKS=50

# Print job information
echo "========================================"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID of $NUM_TASKS"
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
echo "Python version: $(python --version)"
echo ""

# Create chunks directory if it doesn't exist
mkdir -p ../chunks


echo "Input files verified."
echo "Starting task $SLURM_ARRAY_TASK_ID..."
echo "========================================"

# Run the parallel Python script
python errorcorrection_parallel.py \
    --task-id $SLURM_ARRAY_TASK_ID \
    --num-tasks $NUM_TASKS \
    --output-dir ../chunks

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully!"
else
    echo ""
    echo "========================================"
    echo "ERROR: Task $SLURM_ARRAY_TASK_ID failed with exit code $?"
fi

echo "End time: $(date)"
echo "========================================"
