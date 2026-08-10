#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH -n 12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J bcl2fastq
#SBATCH --output=bcl2fastq-%A_%a.out
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# Load Modules -----
module load bcl2fastq2/2.20

# Set working directory
cd /home/groups/hpcbio/temp_fastqs/RNA-Seq/rhodes/2025Nov-scRNASeq/raw-seq/unzip/

# Run bcl2fastq -----
bcl2fastq -R /home/groups/hpcbio/temp_fastqs/RNA-Seq/rhodes/2025Nov-scRNASeq/raw-seq/unzip \
--sample-sheet /home/groups/hpcbio/temp_fastqs/RNA-Seq/rhodes/2025Nov-scRNASeq/raw-seq/unzip/SampleSheet_Rhodes.csv \
--output-dir /home/groups/hpcbio/temp_fastqs/RNA-Seq/rhodes/2025Nov-scRNASeq/raw-seq/Unaligned_Rhodes \
--loading-threads $SLURM_NTASKS \
--processing-threads $SLURM_NTASKS \
--writing-threads $SLURM_NTASKS \
--barcode-mismatches 1 \
--use-bases-mask Y28,I8,Y98 \
--ignore-missing-bcl \
--ignore-missing-filter \
--ignore-missing-positions

# Add this to combine lanes 1 and 2
#--no-lane-splitting
