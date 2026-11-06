#make labelled reference for singleR

#https://portal.brain-map.org/atlases-and-data/rnaseq/human-v1-acc-smart-seq
#smart-seq data

library(SingleCellExperiment)
library(data.table)
library(Matrix)
library(Seurat)
library(scMerge)


#####################
## immune reference
#####################

#make reference using Allen brain institute M1 sample, 
#immune cells from GBM https://www.brainimmuneatlas.org/download.php
#to be run on cluster, needs ~100G
setwd("/data/users/alyne/gliomult/svn/analyse/results")

#read in counts, log normalise (simpler than SCT since there are multiple samples)
#don't think single cell reference needs to be log normalised, but may make more sense when aggregating


#immune cells from GBM recurrence
data_dir = "10x_GBM_immune_rec/filtered_feature_bc_matrix_HumanRecurrentGBM/filtered_gene_bc_matrices/GRCh38"
expr.rec = Read10X(data_dir)
gbm.imm <- CreateSeuratObject(counts = expr.rec, min.cells = 3, min.features = 200)
gbm.imm <- NormalizeData(gbm.imm, normalization.method = "LogNormalize", scale.factor = 10000)

sce2 <- SingleCellExperiment(assays = list(counts = gbm.imm[["RNA"]]$data))
cell.type.labs = fread("10x_GBM_immune_rec/annot_Human_R_GBM_Full.csv.gz",
                       sep = ",", header = T, stringsAsFactors = F, data.table = F)
sce2$label = as.character(sapply(colnames(sce2), function(x) {
  if(x %in% cell.type.labs$cell){
    return(cell.type.labs$cluster[which(cell.type.labs$cell==x)])
  }else{
    return("")
  }
}))
#remove cells with no label
sce2 = sce2[,which(sce2$label!="")]


#healthy brain cells from Allen brain atlas M1 sample
#https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x
counts = fread("allen_brain_map_M1/M1_matrix.csv.gz", data.table = F)
rownames(counts) = counts[,1]
counts = counts[,-1]
health.br <- CreateSeuratObject(counts = t(counts), min.cells = 3, min.features = 200)
health.br <- NormalizeData(health.br, normalization.method = "LogNormalize", scale.factor = 10000)

sce3 <- SingleCellExperiment(assays = list(counts = health.br[["RNA"]]$data))

#read in labels
cell.type.labs = fread("allen_brain_map_M1/M1_metadata.csv.gz",
                       sep = ",", header = T, stringsAsFactors = F, data.table = F)
sce3$label = as.character(sapply(colnames(sce3), function(x) {
  if(x %in% cell.type.labs$sample_name){
    return(cell.type.labs$subclass_label[which(cell.type.labs$sample_name==x)])
  }else{
    return("")
  }
}))


#combine 2 to get full list
sce_list = list(sce2, sce3)
sce.comb = sce_cbind(sce_list,exprs = c("counts"))
sce.comb$label = c(sce2$label,sce3$label)
saveRDS(sce.comb, file = "SingleR_ref_M1_immune.rds")



