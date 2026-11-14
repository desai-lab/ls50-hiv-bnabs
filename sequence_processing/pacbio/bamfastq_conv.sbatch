#!/bin/bash
#SBATCH -p shared
#SBATCH -t 1-0:00
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o bam2fastq_%j.out
#SBATCH -e bam2fastq_%j.err

module load python/3.12.5-fasrc01
mamba activate pbtk_env

bam2fastq -o pacbio_reads /n/ngsdata/pacbio/260123_REV0057/r84147_20260123_183904/2_A01/hifi_reads/m84147_260125_044231_s4.hifi_reads.bam -u