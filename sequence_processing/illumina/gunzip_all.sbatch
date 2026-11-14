#!/bin/bash
#SBATCH --job-name=gunzip_all
#SBATCH --output=../logs/gunzip_%j.out
#SBATCH --error=../logs/gunzip_%j.err
#SBATCH --partition=desai
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00

# Directory containing .gz files (change as needed)
TARGET_DIR="/n/netscratch/desai_lab/Lab/joao/ls50/merged_fastq"

echo "Decompressing all .gz files in: ${TARGET_DIR}"
echo "Start time: $(date)"

cd "${TARGET_DIR}" || { echo "Cannot cd to ${TARGET_DIR}"; exit 1; }

count=0
for f in *.gz; do
    [ -e "$f" ] || { echo "No .gz files found."; exit 0; }
    echo "Decompressing: $f"
    gunzip "$f"
    count=$((count + 1))
done

echo "Done. Decompressed ${count} files."
echo "End time: $(date)"
