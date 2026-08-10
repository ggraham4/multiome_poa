#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=50G
#SBATCH -p hpcbioamd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J preprocess
#SBATCH --output=preprocess-%A_%a.out
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out
#SBATCH --array=1-4%2

# The -D option above it to put the slurm out files in their own directory.

# Load modules ------
module load R/4.4.0-IGB-gcc-8.2.0

# If you need to install packages, you should be on a worker node!!!!

# Run R script ------
# Loop over experiments
echo "Running Sample ${SLURM_ARRAY_TASK_ID}"
Rscript ../Individual_preprocess.R ${SLURM_ARRAY_TASK_ID}
