#!/bin/bash
#SBATCH -J grab_barcodes
#SBATCH -p desai
#SBATCH --cpus-per-task=1
#SBATCH -t 0-4:00
#SBATCH --mem=16G
#SBATCH -o ../logs/grab_barcodes/grabbarcodes_%A_%a.out
#SBATCH -e ../logs/grab_barcodes/grabbarcodes_%A_%a.err
#SBATCH --array=1-188        # Replace N with number of samples

# ============================================
# CUSTOMIZE THESE PATHS - EDIT THESE LINES!
# ============================================
INPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq/"
OUTPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/grab_barcodes/"
SAMPLE_LIST="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq/sample_list2.txt"
# ============================================



module load python
source activate barcode_env

# Create output and log directories
mkdir -p ${OUTPUT_DIR}
mkdir -p '../logs/grab_barcodes'

# Get sample name from array task ID
SAMPLE_FULL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
SAMPLE=$(basename ${SAMPLE_FULL})


echo "========================================"
echo "Sample: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "========================================"

python grab_barcodes.py ${SAMPLE} ${INPUT_DIR} ${OUTPUT_DIR}

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
    echo "ERROR: Grab Barcodes failed for ${SAMPLE}"
    echo "Exit code: ${EXIT_CODE}"
    echo "========================================"
    exit ${EXIT_CODE}
fi

