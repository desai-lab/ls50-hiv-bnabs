#!/bin/bash
#SBATCH -J ec_barcodes_sc
#SBATCH -p desai
#SBATCH --cpus-per-task=8
#SBATCH -t 0-48:00
#SBATCH --mem=32G
#SBATCH -o ../logs/ec_barcodes/rerun5_ecbarcodes_%A_%a.out
#SBATCH -e ../logs/ec_barcodes/rerun5_ecbarcodes_%A_%a.err
#SBATCH --array=1-188

# ============================================
# CUSTOMIZE THESE PATHS - EDIT THESE LINES!
# ============================================
INPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_barcodes/"
OUTPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/umi_counts/"
SAMPLE_LIST="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq/sample_list2.txt"
PACBIO_DIR="/n/holylabs/desai_lab/Lab/joao/ls50/pacbio/"
# ============================================

# Starcode binary location (exported so the python script can find it)
export STARCODE="${HOME}/software/starcode/starcode"

module load python
module load gcc/12.2.0-fasrc01
source activate barcode_env

# Create output and log directories
mkdir -p ${OUTPUT_DIR}
mkdir -p '../logs/ec_barcodes'

# Get sample name from array task ID
SAMPLE_FULL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
SAMPLE=$(basename ${SAMPLE_FULL})

echo "========================================"
echo "Sample: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Output directory: ${OUTPUT_DIR}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "========================================"

# Verify starcode is available
if [ ! -f "${STARCODE}" ]; then
    echo "ERROR: starcode not found at ${STARCODE}"
    exit 1
fi

python ec_count_barcodes_starcode.py ${SAMPLE} ${INPUT_DIR} ${OUTPUT_DIR} ${PACBIO_DIR}

# Check exit status
EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "SUCCESS: ${SAMPLE}"
    echo "Finished: $(date)"
    echo "========================================"
else
    echo ""
    echo "========================================"
    echo "ERROR: ec_count_barcodes failed for ${SAMPLE}"
    echo "Exit code: ${EXIT_CODE}"
    echo "========================================"
    exit ${EXIT_CODE}
fi

