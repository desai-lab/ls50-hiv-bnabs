#!/bin/bash
#SBATCH -J grab_barcodes_undet
#SBATCH -p desai
#SBATCH --cpus-per-task=1
#SBATCH -t 0-24:00
#SBATCH --mem=64G
#SBATCH -o ../logs/grab_barcodes/undet1_%j.out
#SBATCH -e ../logs/grab_barcodes/undet1_%j.err
#SBATCH --mail-type=END,FAIL          # Email notifications on job end or failure
#SBATCH --mail-user=joaoascensao@g.harvard.edu

# ============================================
# CUSTOMIZE THESE PATHS - EDIT THESE LINES!
# ============================================
INPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq/"
OUTPUT_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/grab_barcodes/"

ORIG_FASTQ="/n/ngsdata/20260130_LH00541_0116_A23JVJCLT4_SUB17671/fastq/Undetermined_S0_L003_R1_001.fastq.gz"
PLATE="1"
SAMPLE="Undetermined_S0_L003"

# ============================================

module load python
source activate barcode_env

# Create output and log directories
mkdir -p ${OUTPUT_DIR}
mkdir -p '../logs/grab_barcodes'


echo "========================================"
echo "Sample: ${SAMPLE}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Output prefix: ${OUTPUT_DIR}/${SAMPLE}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "========================================"

python grab_barcodes_undetermined.py ${SAMPLE} ${INPUT_DIR} ${OUTPUT_DIR} ${ORIG_FASTQ} ${PLATE}

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

