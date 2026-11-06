#project each sample onto atlas data
rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
script.dir = as.character(args[3])

library(data.table)
library(dplyr)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(Signac)
library(future)
options(future.globals.maxSize = 100 * 1024 ^ 3)

pt.sz = 0.02

setwd("/data/tmp/fetal_atlas")

#load atlas data and run Seurat pipeline
counts1 = fread("exprMatrix_fetal_filt.tsv.gz", data.table = F, sep="\t")
rownames(counts1) = counts1[,1]
counts1 = counts1[,-1]

counts2 = fread("exprMatrix_adult_filt.tsv.gz", data.table = F, sep="\t")
rownames(counts2) = counts2[,1]
counts2 = counts2[,-1]

#find intersect of rownames
int.genes = sort(intersect(rownames(counts1),rownames(counts2)))
counts1 = counts1[which(rownames(counts1) %in% int.genes),]
counts1 = counts1[int.genes, ]

counts2 = counts2[which(rownames(counts2) %in% int.genes),]
counts2 = counts2[int.genes, ]


counts = cbind(counts1,counts2)

#read in meta data
md1 = read.table("meta_fetal_filt.tsv", header=T, sep="\t", row.names = 1)
md2 = read.table("meta_adult_filt.tsv", header=T, sep="\t", row.names = 1)

#collapse adult neurons into GABAergic and Glutamatergic
md2$Class[which(md2$Class %in% c("LAMP5","VIP","SST","PVALB","PAX6"))] = "GABAergic"
md2$Class[which(md2$Class %in% c("IT","L4 IT","L5 ET","L5/6 IT Car3",
                                 "L5/6 NP","L6 CT","L6b"))] = "Glutamatergic"

md = data.frame(atlas.cell.type = c(paste0(md1$Type.v2,"_fetal"),paste0(md2$Class,"_adult")), 
                Gestational_week = c(md1$Gestational_week,rep("adult",nrow(md2))),
                Dataset = c(md1$Dataset,md2$Dataset))
rownames(md) = c(rownames(md1),rownames(md2))

gest.week.num = sort(as.numeric(unique(md$Gestational_week)[which(!(unique(md$Gestational_week) %in% c("adult","postnatal")))]))
md$Gestational_week = factor(md$Gestational_week, 
                             levels = c(gest.week.num,"postnatal","adult"), 
                             ordered=T)

atlas.fet <- CreateSeuratObject(counts = counts, project = "joint_atlas", min.cells = 10, 
                                min.features = 500, meta.data = md)

atlas.fet[["RNA"]] <- split(atlas.fet[["RNA"]], f = atlas.fet$Dataset)

atlas.fet <- NormalizeData(atlas.fet, normalization.method = "LogNormalize", scale.factor = 10000)
atlas.fet <- FindVariableFeatures(atlas.fet, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(atlas.fet)
atlas.fet <- ScaleData(atlas.fet, features = all.genes)
atlas.fet <- RunPCA(atlas.fet, features = VariableFeatures(object = atlas.fet))
atlas.fet <- FindNeighbors(atlas.fet, dims = 1:30)
atlas.fet <- FindClusters(atlas.fet, resolution = 0.5, cluster.name = "unintegrated_clusters")
atlas.fet <- RunUMAP(atlas.fet, dims = 1:30, reduction.name = "umap.unintegrated")

atlas.fet <- IntegrateLayers(object = atlas.fet, method = RPCAIntegration,
                             orig.reduction = "pca", new.reduction = "integrated.rpca", 
                             verbose = FALSE, k.weight = 20)
atlas.fet <- FindNeighbors(atlas.fet, reduction = "integrated.rpca", dims = 1:30)
atlas.fet <- FindClusters(atlas.fet, resolution = 0.5, cluster.name = "rpca_clusters")
atlas.fet <- RunUMAP(atlas.fet, reduction = "integrated.rpca", dims = 1:30, 
                     reduction.name = "umap.rpca", return.model=T)

#umap plot of atlas
p1 = DimPlot(atlas.fet, group.by = "Gestational_week",
             reduction = 'umap.rpca', pt.size = pt.sz, shuffle=T)
p2 = DimPlot(atlas.fet, group.by = "atlas.cell.type",
             reduction = 'umap.rpca', pt.size = pt.sz, shuffle=T)
print(p1 + p2)
ggsave(paste0(data.dir,"/figures_malig/umap_joint_atlas_Bhaduri_rpca.pdf"), width=17, height=5)

#save atlas
saveRDS(atlas.fet, file=paste0(data.dir,"/joint_atlas_filt.rds"))

print(table(atlas.fet$Gestational_week))

print(table(atlas.fet$atlas.cell.type))

print(paste0("Total number of cells in atlas is ",ncol(atlas.fet)))



