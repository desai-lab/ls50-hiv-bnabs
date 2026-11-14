#!/bin/bash
#SBATCH -J ngmerge_array
#SBATCH -p desai
#SBATCH -c 4
#SBATCH -t 0-2:00
#SBATCH --mem=8G
#SBATCH -o ../logs/ngmerge_%A_%a.out
#SBATCH -e ../logs/ngmerge_%A_%a.err
#SBATCH --array=1-190        # Replace N with number of samples, limit 10 concurrent jobs

# ============================================
# CUSTOMIZE THESE PATHS - EDIT THESE LINES!
# ============================================
OUTPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq"
# ============================================

SAMPLE_LIST="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq/sample_list.txt"

# Create output and log directories
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get sample name from array task ID
SAMPLE_FULL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
SAMPLE=$(basename ${SAMPLE_FULL})

# Construct R1 and R2 filenames
R1="${SAMPLE_FULL}_R1_001.fastq.gz"
R2="${SAMPLE_FULL}_R2_001.fastq.gz"

# Error checking
if [ ! -f "$R1" ]; then
    echo "ERROR: R1 file not found: ${R1}"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "ERROR: R2 file not found: ${R2}"
    exit 1
fi

echo "========================================"
echo "Sample: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Output prefix: ${OUTPUT_DIR}/${SAMPLE}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "========================================"

# Run NGmerge
# -1: Forward reads (R1)
# -2: Reverse reads (R2)
# -o: Output prefix
# -m: Minimum overlap (80bp as specified)
# -n: Number of threads (matches SLURM -c)
# -e: Expected read length (150bp for NovaSeq 2x150)
# -v: Verbose output
# -d: Log dovetailing
# -z: Gzip output files (recommended since inputs are gzipped)

NGmerge -1 ${R1} \
        -2 ${R2} \
        -o ${OUTPUT_DIR}/${SAMPLE} \
        -m 80 \
        -n 4 \
        -e 150 \
        -v \
        -d \
        -z

# Check exit status
EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "SUCCESS: ${SAMPLE}"
    echo "Finished: $(date)"
    echo ""
    echo "Output files created:"
    ls -lh ${OUTPUT_DIR}/${SAMPLE}*
    echo "========================================"
else
    echo ""
    echo "========================================"
    echo "ERROR: NGmerge failed for ${SAMPLE}"
    echo "Exit code: ${EXIT_CODE}"
    echo "========================================"
    exit ${EXIT_CODE}
fi
