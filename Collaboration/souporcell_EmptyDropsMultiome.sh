#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=drnevich@illinois.edu
#SBATCH -J souporcell
#SBATCH --array=1-6
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2024Feb-scRNASeq/src/slurm-out

# The -D option above it to put the slurm out files in their own directory.

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load souporcell/2.5

# List of pools and # of clusters expected
SAMPLEID=`head -n $SLURM_ARRAY_TASK_ID ../pools_clusters.txt | tail -n 1 | cut -f1`
CLUSTERS=`head -n $SLURM_ARRAY_TASK_ID ../pools_clusters.txt | tail -n 1 | cut -f2`

### specifying the output folder
cd /home/groups/hpcbio/RNA-Seq/projects/rhodes/2024Feb-scRNASeq/

### Run souporcell pipeline

souporcell souporcell_pipeline.py \
 -i results/cellranger/s${SAMPLEID}/s${SAMPLEID}_gex_possorted_bam.bam \
 -b results/EmptyDropsMultiome/Pool${SAMPLEID}_Clownfish_eDM.barcodes_PASSED.tsv \
 -t ${SLURM_NTASKS} \
 -f data/references/GCF_022539595.1_ASM2253959v1_genomic.fna \
 -o results/EmptyDropsMultiome/s${SAMPLEID} \
 -k $CLUSTERS
 
