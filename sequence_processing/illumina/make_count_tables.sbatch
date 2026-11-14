#!/bin/bash
#SBATCH -J make_count_tables
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00
#SBATCH --mem=8G
#SBATCH -o ../logs/make_count_tables_%j.out
#SBATCH -e ../logs/make_count_tables_%j.err

module load python
source activate barcode_env

python make_count_tables.py

