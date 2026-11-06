#!/bin/bash
#SBATCH --job-name=diff_peak_ini_rec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200gb
#SBATCH --time=10:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr


#script to launch DAR analysis for paired ini-rec oligodendroglioma

source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/archr_env
# module load R

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6

SCRIPT_DIR="/mnt/beegfs/home/alyne/persistent/ArchR_scripts"

cancer="oligo"

$CONDA_PREFIX/bin/Rscript "$SCRIPT_DIR"/diff_analysis_peaks_LIBRA.r \
  --archr_path /mnt/beegfs/home/alyne/stageout/ArchR_output_"$cancer"/project2 \
  --out_directory /mnt/beegfs/home/alyne/stageout/ArchR_output_"$cancer"/diff_output_merged \
  --diff_annotation ini_rec \
  --group_forward rec1 \
  --group_background ini \
  --cell_type_column cell_type \
  --replicate_column Patient.number \
  --sample_column Sample.number \
  --paired TRUE \
  --de_method DESeq2 \
  --de_type Wald \
  --binarization FALSE \
  --min_cell_threshold 50 \
  --max_cells 5000 \
  --reproducibility "max(ceiling(n/3),2)" \
  --exclude_patterns "rec2" \
  --seed 123 \
  --motifdb JASPAR2024 \
  --threshold_FDR 0.05 \
  --threshold_Log2FC 0.5 \
  --peak_call FALSE


