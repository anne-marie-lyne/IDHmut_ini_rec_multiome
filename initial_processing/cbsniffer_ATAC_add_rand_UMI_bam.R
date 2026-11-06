#script to make figures summarising whole cohort

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
samp = as.character(args[1])
tmp.dir = as.character(args[2])
nReads = as.numeric(args[3])

library(stringi)
library(data.table)


#can't figure out how to run cb_sniffer on ATAC data without adding "fake" UMI to bam file
#this script produces text file with required number of UMIs to be added to ATAC bam file
#output is in correct format for bam, similar to gex bam file from cellranger

#simulate random UMI barcodes to add to bam file
bc.vec = stringi::stri_rand_strings(n = nReads, length = 12, pattern = "[ATCG]")

#format to be like 10x barcodes
output.vec = paste0("UB:Z:",bc.vec)

fwrite(list(output.vec), file = paste0(tmp.dir,"/umis.txt"), row.names = F, col.names = F, quote=F)















