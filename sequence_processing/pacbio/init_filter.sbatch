#!/bin/bash
#SBATCH -p desai
#SBATCH -t 10-0:00
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o ../logs/initfilter_%j.out
#SBATCH -e ../logs/initfilter_%j.err

bash init_filter.sh ../pacbio_reads.fastq
