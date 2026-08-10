#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=140G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J CR-count
#SBATCH --output=CR-count-%A_%a.out
#SBATCH --array=1-4%2
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out

# The -D option above it to put the slurm out files in their own directory.
# Now go up one to the src directory, where files.txt is
cd ../

# Pull sample names and fastq locations and put in variables:
Sample_ID=`head files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
Sample_name=`head files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`
File_folder=`head files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/9.0.0

### specifying the output folder
cd ../results/cellranger

### point to a cellranger reference
DB=/home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/data/references/STAR_Ref_Genome_ocellaris

### Run app on file
cellranger count --id=${Sample_ID} \
    --localcores=$SLURM_NTASKS  \
    --localmem 128 \
    --transcriptome=$DB \
    --fastqs=$File_folder \
    --sample=$Sample_name \
    --create-bam=false

# Rename web summary and cloupe files to add sample ID
mv $Sample_ID/outs/web_summary.html $Sample_ID/outs/${Sample_ID}_web_summary.html
mv $Sample_ID/outs/cloupe.cloupe $Sample_ID/outs/${Sample_ID}.cloupe


