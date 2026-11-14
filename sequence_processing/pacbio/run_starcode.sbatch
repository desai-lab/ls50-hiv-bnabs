#!/bin/bash
#SBATCH --job-name=starcode
#SBATCH --partition=desai          # Use 'test' for quick tests (<1hr), 'shared' for longer jobs
#SBATCH --time=0-04:00:00           # 4 hours (adjust based on data size)
#SBATCH --mem=32G                   # 32GB RAM (starcode can be memory-intensive for large datasets)
#SBATCH --cpus-per-task=8           # Number of threads for parallel processing
#SBATCH --output=../logs/starcode_%j.out    # Standard output log (%j = job ID)
#SBATCH --error=../logs/starcode_%j.err     # Standard error log
#SBATCH --mail-type=END,FAIL        # Email notifications (optional)
#SBATCH --mail-user=joaoascensao@g.harvard.edu

#==============================================================================
# STARCODE BARCODE CLUSTERING
# 
# This script runs starcode for DNA barcode error correction/clustering.
# Starcode finds all pairs of sequences within a given Levenshtein distance
# and clusters them together.
#==============================================================================

# Exit on any error
set -e

# Load required modules
module load gcc/12.2.0-fasrc01

#------------------------------------------------------------------------------
# CONFIGURATION - MODIFY THESE VARIABLES
#------------------------------------------------------------------------------

# Path to starcode executable
STARCODE="${HOME}/software/starcode/starcode"

# Input file: one barcode per line (or FASTA/FASTQ)
INPUT_FILE="../barcodes.txt"

# Output file
OUTPUT_FILE="../clustered_barcodes.txt"

# Maximum Levenshtein distance for merging (your max_mismatch_bc parameter)
MAX_DISTANCE=6

# Number of threads (should match --cpus-per-task above)
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Clustering algorithm:
#   "mp"     = Message Passing (default) - iterative bottom-up clustering
#   "sphere" = Sphere clustering - greedy, largest sequences absorb matches
#   "cc"     = Connected Components - groups all sequences within distance d
#              (most similar to your current NetworkX approach)
CLUSTER_ALG="cc"

#------------------------------------------------------------------------------
# VALIDATION
#------------------------------------------------------------------------------

echo "=========================================="
echo "STARCODE CLUSTERING JOB"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Started: $(date)"
echo ""

# Check that starcode exists
if [ ! -f "${STARCODE}" ]; then
    echo "ERROR: Starcode not found at ${STARCODE}"
    echo "Please install starcode first. See installation instructions."
    exit 1
fi

# Check that input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# Count input sequences
INPUT_COUNT=$(wc -l < "${INPUT_FILE}")
echo "Input file: ${INPUT_FILE}"
echo "Number of sequences: ${INPUT_COUNT}"
echo "Max Levenshtein distance: ${MAX_DISTANCE}"
echo "Threads: ${THREADS}"
echo "Clustering algorithm: ${CLUSTER_ALG}"
echo ""

#------------------------------------------------------------------------------
# RUN STARCODE
#------------------------------------------------------------------------------

echo "Starting starcode..."
echo "------------------------------------------"

# Build the starcode command based on clustering algorithm
STARCODE_CMD="${STARCODE} -d ${MAX_DISTANCE} -t ${THREADS} -i ${INPUT_FILE} -o ${OUTPUT_FILE} --print-clusters"

case ${CLUSTER_ALG} in
    "sphere"|"s")
        STARCODE_CMD="${STARCODE_CMD} --spheres"
        ;;
    "cc"|"connected")
        STARCODE_CMD="${STARCODE_CMD} --connected-comp"
        ;;
    "mp"|"message")
        # Message passing is default, no extra flag needed
        ;;
    *)
        echo "WARNING: Unknown clustering algorithm '${CLUSTER_ALG}', using default (message passing)"
        ;;
esac

echo "Command: ${STARCODE_CMD}"
echo ""

# Run starcode and time it
START_TIME=$(date +%s)

${STARCODE_CMD}

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "------------------------------------------"
echo "Starcode completed in ${ELAPSED} seconds"

#------------------------------------------------------------------------------
# SUMMARIZE RESULTS
#------------------------------------------------------------------------------

if [ -f "${OUTPUT_FILE}" ]; then
    CLUSTER_COUNT=$(wc -l < "${OUTPUT_FILE}")
    echo ""
    echo "=========================================="
    echo "RESULTS SUMMARY"
    echo "=========================================="
    echo "Input sequences: ${INPUT_COUNT}"
    echo "Output clusters: ${CLUSTER_COUNT}"
    echo "Reduction: $(echo "scale=2; (1 - ${CLUSTER_COUNT}/${INPUT_COUNT}) * 100" | bc)%"
    echo ""
    echo "Output file: ${OUTPUT_FILE}"
    echo ""
    echo "First 10 clusters:"
    head -10 "${OUTPUT_FILE}"
    echo ""
else
    echo "ERROR: Output file was not created"
    exit 1
fi

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="

