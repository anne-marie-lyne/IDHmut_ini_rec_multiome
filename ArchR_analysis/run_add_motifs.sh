#!/bin/bash
#SBATCH --job-name=diff_peak_PCV_oligo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=300gb
#SBATCH --time=48:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr

source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/archr_env
# module load R

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6

SCRIPT_DIR="/mnt/beegfs/home/alyne/persistent/ArchR_scripts"

cancer="oligo"

$CONDA_PREFIX/bin/Rscript "$SCRIPT_DIR"/add_motif.R \
