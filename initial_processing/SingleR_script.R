#Anne-Marie Lyne 01/02/2021

#script to run SingleR on Neftel et al (2019) 10x dataset

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
wd = as.character(args[1])
pbs_arrayid = as.numeric(args[2])
ref.path.comb = as.character(args[3])
samp.file = as.character(args[4])
ref = as.character(args[5])


library(SingleR)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

setwd(wd)

#get list of all samples
sample.info = read.table(samp.file, header = F, stringsAsFactors = F, fill = T,
                         colClasses = c("character","character","numeric"))
samples = sample.info[,1]
samples.newnames = sample.info[,2]


#load Seurat object
oligo = readRDS(paste0("../oligo_filt_", samples.newnames[pbs_arrayid], ".rds"))
samp = unique(oligo$orig.ident)


#load reference (made using script SingleR_make_SCE_reference.R)
#NB ref has been log normalised
#not necessary for single cell ref when wilcox used,
#but may be necessary when single cell ref is aggregated by cell type
sce.comb = readRDS(ref.path.comb)


pred <- SingleR(test = oligo[["RNA"]]$counts, ref = counts(sce.comb),
                labels = sce.comb$label,
                de.method = "wilcox",
                aggr.ref = T,
                fine.tune = T,
                prune = T)

#write predictions to file so can be added to seurat object afterwards
write.table(pred, file = paste0("SingleR_",ref,"_",samp,".csv"),
       quote = F, row.names = T, col.names = T, sep = ",")


#make heatmap of scores for each cell for all cell types
plotScoreHeatmap(pred, filename = paste0("../../figures/SingleR_score_hm_",ref,"_",samp,".pdf"))

#summary of number of cells with label NA
summary(is.na(pred$pruned.labels))

#plot of delta distributions (i.e. difference between top score and second best)
plotDeltaDistribution(pred, ncol = 3)
ggsave(paste0("../../figures/SingleR_delta_",ref,"_",samp,".pdf"))


