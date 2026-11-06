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

#make short version of variant file for test
# head -n 2 "$DATA_DIR"/TCGA_LGG_variants.txt > "$DATA_DIR"/TCGA_LGG_variants_short.txt

cd "$DATA_DIR"/analysis_"$SAMPLE"

mkdir cb_sniffer_"$SAMPLE"
cd cb_sniffer_"$SAMPLE"
# rm *_ATAC*
# 
# 
# 
#find location on KDI of atac and rna bams
ATAC_LOC=$(echo "$KDI_ATAC" | awk -v l="$BATCH" 'BEGIN {FS=";"} {print $l}')
ATAC_BAM="$ATAC_LOC"/"$SAMPLE_ORIG"/bam

#move bam and index to /data/tmp
mkdir /data/tmp/cb_snif_tmp
mkdir /data/tmp/cb_snif_tmp/"$SAMPLE"
cp "$ATAC_BAM"/atac_possorted_bam.bam /data/tmp/cb_snif_tmp/"$SAMPLE"/
cp "$ATAC_BAM"/atac_possorted_bam.bam.bai /data/tmp/cb_snif_tmp/"$SAMPLE"/

echo "ATAC files copied"

cd /data/tmp/cb_snif_tmp/"$SAMPLE"


# #count number of reads in bam file
NREAD=$(samtools view atac_possorted_bam.bam | wc -l)

echo "reads counted"

#run R code to simulate this number of random UMIs
R_FILE="${SCRIPT_DIR}/cbsniffer_ATAC_add_rand_UMI_bam.R"
$R_EXEC $R_FILE "$SAMPLE" "/data/tmp/cb_snif_tmp/${SAMPLE}" "$NREAD"

echo "random UMIs simulated"


#write header to new file
samtools view -H atac_possorted_bam.bam > atac_possorted_bam_umi.txt
#join umis to each row in bam file
paste <(samtools view atac_possorted_bam.bam) umis.txt >> atac_possorted_bam_umi.txt
#gzip and then convert to bam
gzip atac_possorted_bam_umi.txt
gzip -dc atac_possorted_bam_umi.txt.gz | samtools view -Sb -o atac_possorted_bam_umi.bam -

echo "new bam file made"


#index
cd "$DATA_DIR"/analysis_"$SAMPLE"/cb_sniffer_"$SAMPLE"
samtools index /data/tmp/cb_snif_tmp/"$SAMPLE"/atac_possorted_bam_umi.bam


#run cb_sniffer
gzip -d "$DATA_DIR"/analysis_"$SAMPLE"/cells_pass_filt_both_"$SAMPLE".txt.gz
python3.6 "$SCRIPT_DIR"/cb_sniffer.py /data/tmp/cb_snif_tmp/"$SAMPLE"/atac_possorted_bam_umi.bam "$DATA_DIR"/TCGA_LGG_variants.txt "$DATA_DIR"/analysis_"$SAMPLE"/cells_pass_filt_both_"$SAMPLE".txt "$SAMPLE"_ATAC


gzip "$DATA_DIR"/analysis_"$SAMPLE"/cells_pass_filt_both_"$SAMPLE".txt








