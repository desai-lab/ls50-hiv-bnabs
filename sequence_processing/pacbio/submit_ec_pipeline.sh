#!/bin/bash
# ============================================================================
# CONVENIENCE SCRIPT: Submit parallel error correction pipeline
# ============================================================================
# This script submits the array job and automatically schedules the merge
# job to run after all array tasks complete.
#
# Usage: ./submit_pipeline.sh
# ============================================================================

echo "Submitting parallel barcode error correction pipeline..."
echo ""

# Create necessary directories
mkdir -p ../chunks
mkdir -p ../logs

# Submit array job and capture job ID
ARRAY_JOB_ID=$(sbatch --parsable errorcorrection_array.sbatch)

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to submit array job!"
    exit 1
fi

echo "Submitted array job: $ARRAY_JOB_ID"

# Submit merge job with dependency on array job completion
MERGE_JOB_ID=$(sbatch --parsable --dependency=afterok:$ARRAY_JOB_ID errorcorrection_merge.sbatch)

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to submit merge job!"
    echo "Array job $ARRAY_JOB_ID is still running. Submit merge manually after completion."
    exit 1
fi

echo "Submitted merge job: $MERGE_JOB_ID (will run after array job completes)"
echo ""
echo "========================================"
echo "Pipeline submitted successfully!"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo ""
echo "View array job logs:"
echo "  tail -f ../logs/barcode_ec_${ARRAY_JOB_ID}_*.out"
echo ""
echo "View merge job logs (after array completes):"
echo "  tail -f ../logs/barcode_ec_merge_${MERGE_JOB_ID}.out"
echo "========================================"
