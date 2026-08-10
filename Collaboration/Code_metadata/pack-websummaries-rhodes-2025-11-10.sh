#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J pack-html
#SBATCH --output=pack-html-%A.out
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out

# pick cellranger or spaceranger
cd ../../results/cellranger
# cd ../../results/spaceranger

# Find all web summary files and zip

find */outs -type f -name "*web_summary.html" | zip -j web_summaries.zip -@ 


