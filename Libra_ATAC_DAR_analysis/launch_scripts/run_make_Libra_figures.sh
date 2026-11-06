#!/bin/bash
#SBATCH --job-name=make_figs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=2:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr

#script to launch code to make figures
#can also just directly use R script analyse_Libra_output.R
#and set parameters in script

source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/archr_env
# module load R

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6

SCRIPT_DIR="/mnt/beegfs/home/alyne/persistent/ArchR_scripts"

cancer="oligo"
folder="diff_output_merged"
str_ana="grade"
gr_for="3"
gr_back="2"
analysis_name="oligodendroglioma grade"
archr_proj="project2"

$CONDA_PREFIX/bin/Rscript "$SCRIPT_DIR"/analyse_Libra_output.R $folder $str_ana $cancer "$analysis_name" $archr_proj $gr_for $gr_back


# folder="diff_output"
# archr_proj="project1"
# 
# $CONDA_PREFIX/bin/Rscript "$SCRIPT_DIR"/analyse_Libra_output.R $folder $str_ana $cancer "$analysis_name" $archr_proj $gr_for $gr_back






