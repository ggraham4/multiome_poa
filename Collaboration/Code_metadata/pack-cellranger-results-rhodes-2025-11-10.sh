#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -p hpcbioamd,hpcbio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=valizad2@illinois.edu
#SBATCH -J pack-CR
#SBATCH --output=pack-CR-%A_%a.out
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/src/slurm-out
#SBATCH --array=1-5%2

#--array should be equal to the number of directories to be packed either:
# 1. number of samples + Aggregated
# or
# 2. number of pools + Aggregated

# For cellranger multi runs with multiplexed samples, each sample will get a different 
# .zip file while there will be one .zip file for the multiplexed pool "raw" results

# codes below written with the help of ChatGPT

set -Eeuo pipefail
IFS=$'\n\t'

# ===== Dry-run switch =====
# DRYRUN=1 → only print commands
# DRYRUN=0 → actually run commands
DRYRUN=0
run() {
  if [ "${DRYRUN}" -eq 1 ]; then
    echo "[DRYRUN] $*"
  else
    eval "$@"
  fi
}

cd ../../results/cellranger

# 1) Pull directory names by number (or use option 2 below)
i=$(printf "%s\n" */ | sed 's:/$::' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Guard against empty/bad $i (e.g., array index > number of dirs)
if [ -z "${i:-}" ] || [ ! -d "$i" ]; then
  echo "No directory for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} — skipping."
  exit 0
fi

# 2) Alternative: pick from an explicit list
# i=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" packdirnames.txt)

# ======================
# Branch by output type
# ======================

# 1) cellranger multi: detect by per_sample_outs directory FIRST
if [ -d "$i/outs/per_sample_outs" ]; then
  for sampledir in "$i/outs/per_sample_outs"/*; do
    [ -d "$sampledir" ] || continue
    samplename=$(basename "$sampledir")

    # Rename generic web_summary.html -> ${samplename}_web_summary.html
    while IFS= read -r -d '' f; do
      dir=$(dirname "$f")
      if [ -f "$dir/web_summary.html" ]; then
        run mv -- "$dir/web_summary.html" "$dir/${samplename}_web_summary.html"
      fi
    done < <(find "$sampledir" -type f -name 'web_summary.html' -print0)

    # Rename sample_cloupe.cloupe -> ${samplename}.cloupe
    while IFS= read -r -d '' c; do
      d=$(dirname "$c")
      if [ -f "$c" ]; then
        run mv -- "$c" "$d/${samplename}.cloupe"
      fi
    done < <(find "$sampledir" -type f -name 'sample_cloupe.cloupe' -print0)

    # Per-sample zip at top level, exclude BAMs
    run zip -qr "./${i}_${samplename}.zip" "$sampledir/" -x '*bam*'

    # Per-pool zip of multi raw outputs
    run zip -qr "./${i}_multi.zip" "./${i}/outs/multi/" -x '*bam*'
    
  done

# 2) cellranger count/aggr: any *web_summary.html under outs/
elif find "$i/outs" -type f -name '*web_summary.html' -print -quit >/dev/null 2>&1; then
  # Rename every generic web_summary.html -> "${i}_web_summary.html"
  while IFS= read -r -d '' f; do
    dir=$(dirname "$f")
    if [ -f "$dir/web_summary.html" ]; then
      run mv -- "$dir/web_summary.html" "$dir/${i}_web_summary.html"
    fi
  done < <(find "$i/outs" -type f -name 'web_summary.html' -print0)

  # Rename any cloupe.cloupe -> "${i}.cloupe"
  while IFS= read -r -d '' c; do
    d=$(dirname "$c")
    run mv -- "$c" "$d/${i}.cloupe"
  done < <(find "$i" -type f -name 'cloupe.cloupe' -print0)

  # Zip outs/ if present; otherwise zip the whole sample dir (quiet), exclude BAMs
  if [ -d "./$i/outs" ]; then
    run zip -qr "./${i}.zip" "./$i/outs/" -x '*bam*'
  else
    run zip -qr "./${i}.zip" "./$i/" -x '*bam*'
  fi
fi

# ----
# Legacy tar.gz kept for reference (disabled)
# module load pigz
# if test -d "./${i}/outs" ; then
#   tar -cf - --exclude=*bam* ./${i}/outs/ | pigz -p $SLURM_NTASKS > ${i}.tar.gz
# else
#   tar -cf - --exclude=*bam* ./${i}/ | pigz -p $SLURM_NTASKS > ${i}.tar.gz
# fi
