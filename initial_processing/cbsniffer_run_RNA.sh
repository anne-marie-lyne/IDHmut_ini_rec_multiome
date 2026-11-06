#!/bin/bash
#  run_cb_sniffer.sh
#
#  Created by A-M Lyne on 07/12/2021.
#  
# run cb sniffer to detect known variants in 10x data
# 

SAMPLE_ORIG=$(awk -v l="$PBS_ARRAYID" 'NR==l {print $1}' "$SAMPLE_FILE")
SAMPLE=$(awk -v l="$PBS_ARRAYID" 'NR==l {print $2}' "$SAMPLE_FILE")
BATCH=$(awk -v l="$PBS_ARRAYID" 'NR==l {print $4}' "$SAMPLE_FILE")
echo $BATCH
echo "$SAMPLE"


cd "$DATA_DIR"/analysis_"$SAMPLE"

gzip -d cells_pass_filt_both_"$SAMPLE".txt.gz

mkdir cb_sniffer_"$SAMPLE"
cd cb_sniffer_"$SAMPLE"
rm *


#find location on KDI of rna bams
RNA_LOC=$(echo "$KDI_RNA" | awk -v l="$BATCH" 'BEGIN {FS=";"} {print $l}')
RNA_BAM="$RNA_LOC"/"$SAMPLE_ORIG"/bam


source ~/gliomult/svn/analyse/results/miniconda3/etc/profile.d/conda.sh
conda activate cbsniffer_env



if [ -f "$RNA_BAM"/gex_possorted_bam.bam ]; then

python "$SCRIPT_DIR"/cb_sniffer.py "$RNA_BAM"/gex_possorted_bam.bam "$DATA_DIR"/cbsniffer_TCGA_LGG_variants.txt "$DATA_DIR"/analysis_"$SAMPLE"/cells_pass_filt_both_"$SAMPLE".txt "$SAMPLE"_RNA

else

echo "gex bam not in this location, stopping."

fi



cd ..
gzip cells_pass_filt_both_"$SAMPLE".txt










