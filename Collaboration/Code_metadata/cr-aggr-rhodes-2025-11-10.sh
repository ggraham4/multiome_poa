#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=150G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J CR-aggr
#SBATCH --output=CR-aggr-%A.out
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out



### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/9.0.0

### specifying the output folder
cd ../../results/cellranger

cellranger aggr --id=Aggregated \
                  --csv=../../src/AggrList1.csv \
                  --normalize=none \
                  --localcores=$SLURM_NTASKS \
                  --localmem 128

