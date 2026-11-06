rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
pbs.arrayid = as.numeric(args[3])
script.dir = as.character(args[4])

setwd(data.dir)

source(paste0(script.dir,"/make_cell_type_palette.R"))

library(Seurat)
library(Signac)
library(SingleR)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(reticulate)
library(future)
options(future.globals.maxSize = 100 * 1024 ^ 3)


#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric"),
                         stringsAsFactors = F, fill = T)
samples = sample.info[,1]
samples.newnames = sample.info[,2]

iCNV.thresh = 0.99


#############################
## separately run SingleR with SingleR_script.R
## using Allen/imm references
#############################

cl.res = 0.8
# #read in previously saved oligo.filt and add SingleR annotations
oligo.filt = readRDS(paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                            "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))

ref = "imm"

#remove previous metadata columns for this reference, if exists
inds = which(grepl(paste0(".",ref),colnames(oligo.filt@meta.data),fixed=T))
if(length(inds)>0){
  oligo.filt@meta.data = oligo.filt@meta.data[,-inds]
}


#path to SingleR output for this sample
pth = paste0("analysis_",samples.newnames[pbs.arrayid],"/SingleR_",
             samples.newnames[pbs.arrayid],"/SingleR_",ref,"_",samples.newnames[pbs.arrayid],".csv.gz")


#add SingleR output for each cell
if(file.exists(pth)){
  SR = fread(pth, data.table = F)
  cell.type.sr = unlist(sapply(colnames(oligo.filt@assays$SCT),
                               function(y) {SR$pruned.labels[which(SR$V1==y)]}))
  names(cell.type.sr) = names(oligo.filt$orig.ident)
  
  #add to metadata
  oligo.filt@meta.data = cbind(oligo.filt@meta.data,cell.type.sr)
  colnames(oligo.filt@meta.data) = c(colnames(oligo.filt@meta.data)[1:(ncol(oligo.filt@meta.data)-1)], paste0("cell.type.sr.",ref))
  
  #combine cell types
  oligo.filt@meta.data[,ncol(oligo.filt@meta.data)] = combine_cell_types_imm(oligo.filt@meta.data[,ncol(oligo.filt@meta.data)])
  
  #find most common cell type per cluster and percentage
  md.col = which(colnames(oligo.filt@meta.data)==paste0("SCT_snn_res.",cl.res))
  #get unique clusters
  cl.uni = unique(oligo.filt@meta.data[,md.col])
  #find most common annotation for each cluster
  cl.annot = sapply(cl.uni, function(y) {
    annots.tab = table(oligo.filt@meta.data[oligo.filt@meta.data[,md.col]==y,ncol(oligo.filt@meta.data)])
    return(names(annots.tab)[which.max(annots.tab)])
  })
  #find percentage of cells with this annotation
  cl.perc = sapply(cl.uni, function(y) {
    annots.tab = table(oligo.filt@meta.data[oligo.filt@meta.data[,md.col]==y,ncol(oligo.filt@meta.data)])
    return(100*annots.tab[which.max(annots.tab)]/sum(annots.tab))
  })
  #assign most abundant annotation to cluster
  #(if statement there to allow assignment only if x% of cells has same annotation)
  cell.type.sr.cl = sapply(oligo.filt@meta.data[,md.col], function(y) {
    cl.annot[which(cl.uni==y)]
  })
  type.perc.sr.cl = sapply(oligo.filt@meta.data[,md.col], function(y) {
    cl.perc[which(cl.uni==y)]
  })
  names(cell.type.sr.cl) = names(oligo.filt$orig.ident)
  names(type.perc.sr.cl) = names(oligo.filt$orig.ident)
  oligo.filt@meta.data = cbind(oligo.filt@meta.data,cell.type.sr.cl,type.perc.sr.cl)
  colnames(oligo.filt@meta.data) = c(colnames(oligo.filt@meta.data)[1:(ncol(oligo.filt@meta.data)-2)],
                                     paste0("cell.type.sr.cl.",ref),paste0("type.perc.sr.cl.",ref))
  
  # save output separately
  saveRDS(oligo.filt, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                                    "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))
  
  #RNA umap with SingleR ids, SR ids by cluster and CNV info
  p1 = DimPlot(oligo.filt, group.by = paste0("SCT_snn_res.",cl.res), reduction = 'umap.rna', pt.size = 0.005) +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p2 = DimPlot(oligo.filt, group.by = paste0("cell.type.sr.",ref), reduction = 'umap.rna', pt.size = 0.005) +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p3 = FeaturePlot(oligo.filt, features = paste0("type.perc.sr.cl.",ref), reduction = 'umap.rna', pt.size = 0.005) +
    scale_color_gradientn(colors = topo.colors(30), na.value = "gray90") +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p4 = DimPlot(oligo.filt, group.by = paste0("cell.type.sr.cl.",ref), reduction = 'umap.rna', pt.size = 0.005) +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  print(p1 + p2 + p3 + p4 + plot_layout(ncol = 2))
  ggsave(paste0("figures/umaps_",ref,"_celltype_", samples.newnames[pbs.arrayid], ".pdf"), width = 15, height = 8)
}





##############################
## Run code to make merged reference for infercnv (infercnv_make_reference.R)
## Run inferCNV, non-hmm (infercnv_script.R)
## Run cbsniffer
## Then identify malignant clusters based on expression of 1p genes for oligo
## and CNA and SNV for astro using identify_malignant_cells.R
###################################



#############
#### project cells onto adult/fetal atlas then add labels to seurat object
#### use atlas_prelim_processing_fetal_adult.R to process atlas
#############


if(file.exists(paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                      "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))){
  pt.sz = 0.02
  #load atlas
  atlas.fet = readRDS("joint_atlas_filt.rds")
  
  #make new directory
  new.dir = paste0("analysis_",samples.newnames[pbs.arrayid],"/atlas_proj_joint_",samples.newnames[pbs.arrayid])
  dir.create(new.dir)
  
  #read in each sample and project onto atlas
  oligo = readRDS(paste0("analysis_", samples.newnames[pbs.arrayid], "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))
  cancer = unlist(strsplit(oligo$orig.ident[1],"_"))[1]
  DefaultAssay(oligo) = "RNA"
  
  #log normalise to match atlas
  oligo <- NormalizeData(oligo, normalization.method = "LogNormalize", scale.factor = 10000)
  oligo <- FindVariableFeatures(oligo, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(oligo)
  oligo <- ScaleData(oligo, features = all.genes)
  oligo <- RunPCA(oligo, features = VariableFeatures(object = oligo))
  oligo <- FindNeighbors(oligo, dims = 1:30)
  oligo <- FindClusters(oligo, resolution = 0.5)
  oligo <- RunUMAP(oligo, dims = 1:50)
  
  
  #find anchors and transfer labels
  anchors <- FindTransferAnchors(reference = atlas.fet, query = oligo, dims = 1:30,
                                 reference.reduction="integrated.rpca", k.filter = NA)
  oligo <- MapQuery(anchorset = anchors, reference = atlas.fet, query = oligo,
                    refdata = atlas.fet$atlas.cell.type, reference.reduction = "pca", reduction.model = "umap.rpca")
  
  #also transfer gestational week
  predictions <- TransferData(anchorset = anchors, refdata = atlas.fet$Gestational_week, dims = 1:30)
  oligo$pred.id.gest_week = predictions$predicted.id
  
  
  #merge atlas and projected sample to plot joint umaps
  merged.ds = merge(atlas.fet,oligo, merge.data = F, merge.dr = F)
  atlas_proj = factor(c(rep("atlas",ncol(atlas.fet)), rep(oligo$orig.ident[1],ncol(oligo))))
  names(atlas_proj) = colnames(merged.ds)
  merged.ds$atlas_proj = atlas_proj
  merged.ds[["umap_comb"]] <- CreateDimReducObject(embeddings = rbind(Embeddings(atlas.fet, reduction = "umap.rpca")[,1:2],
                                                                      Embeddings(oligo, reduction = "ref.umap")[,1:2]),
                                                   key = "umapcomb_", assay = DefaultAssay(merged.ds))
  
  #plot umaps
  p1 <- DimPlot(atlas.fet, reduction = "umap.rpca", group.by = "Gestational_week") + ggtitle("Atlas - gestational week")
  p2 <- DimPlot(atlas.fet, reduction = "umap.rpca", group.by = "atlas.cell.type", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Atlas - cell type")
  p3 <- DimPlot(merged.ds, reduction = "umap_comb", group.by = "atlas_proj") + ggtitle("Origin")
  p4 <- DimPlot(oligo, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Transferred labels")
  p5 <- DimPlot(oligo, reduction = "ref.umap", group.by = "cell.type.sr.imm") + ggtitle("SingleR")
  p.full = (p1 + p2 + p3)/(p4 + p5)
  ggsave(paste0("figures/umap_joint_atlas_",samples.newnames[pbs.arrayid],".pdf"), plot = p.full, width=13, height=15)
  
  
  #save transferred labels and gestational week
  to.save = data.frame(oligo$predicted.id, oligo$pred.id.gest_week)
  rownames(to.save) = names(oligo$predicted.id)
  fwrite(to.save, file = paste0(new.dir,"/transferred_atlas_labels_",samples.newnames[pbs.arrayid],".tsv.gz"),
         sep="\t", append=F, row.names = T)
  
  
  #switch back to SCT and add labels to Seurat object
  DefaultAssay(oligo) = "SCT"
  atlas.label = sapply(names(oligo$orig.ident), function(x) {to.save$oligo.predicted.id[which(rownames(to.save)==x)]})
  gest_week = sapply(names(oligo$orig.ident), function(x) {to.save$oligo.pred.id.gest_week[which(rownames(to.save)==x)]})
  names(atlas.label) = names(oligo$orig.ident)
  names(gest_week) = names(oligo$orig.ident)
  oligo$atlas.label = atlas.label
  oligo$gest_week = gest_week

  saveRDS(oligo, paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                        "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))
}



##############################
## Then run make_figures_cancer_all_cells.R to merge cells across samples
## and detect remaining doublets and re-annotate cells from less abundant cell types
###################################


#load sample
oligo.filt = readRDS(paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                            "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))

#first remove extra doublets
cancer = unlist(strsplit(oligo.filt$orig.ident[1],"_"))[1]

#read in list of cells to remove and list of cells to annotate
cells.remove = read.table(paste0("probable_doublets_",cancer,".txt"), header=F, row.names = NULL)[,2]
cells.names.rem = sapply(cells.remove, function(x) {
  paste(unlist(strsplit(x, "-"))[1:2],collapse="-")
})
sample.names.rem = sapply(cells.remove, function(x) { unlist(strsplit(x, "-"))[3] })
cells.rem.thissamp = cells.names.rem[which(sample.names.rem==oligo.filt$orig.ident[1])]

#read in cell types after merge
df.cell.types = read.table(paste0("cell_types_after_merge_",cancer,".txt"), header = F)
#get cell names and sample anme separately
df.cell.types.cells = sapply(df.cell.types[,1], function(x) {
  paste(unlist(strsplit(x, "-"))[1:2],collapse="-")
})
df.cell.types.samps = sapply(df.cell.types[,1], function(x) {
  unlist(strsplit(x, "-"))[3]
})
df.cell.types = cbind(df.cell.types, df.cell.types.cells, df.cell.types.samps)
df.cell.types.thissamp = df.cell.types[which(df.cell.types$df.cell.types.samps==oligo.filt$orig.ident[1]),]

#filter to take only non-doublet cells
oligo.filt2 = subset(oligo.filt, cells = names(oligo.filt$orig.ident)[which(!(names(oligo.filt$orig.ident) %in%
                                                                               cells.rem.thissamp))])

#change cell type of cells which need changing
#NB these are generally cell types with very low abundance which need the full dataset to form a cluster
oligo.filt2$cell.type.sr.cl.imm.postmerge = oligo.filt2$cell.type.sr.cl.imm
oligo.filt2$cell.type.sr.cl.imm.postmerge = sapply(seq(length(oligo.filt2$orig.ident)), function(x) {
  if(names(oligo.filt2$orig.ident)[x] %in% df.cell.types.thissamp$df.cell.types.cells){
    return(df.cell.types.thissamp[which(df.cell.types.thissamp$df.cell.types.cells==names(oligo.filt2$orig.ident)[x]),4])
  }else{
    return(oligo.filt2$cell.type.sr.cl.imm.postmerge[x])
  }
})

#ensure these cells are labelled as non-malignant
oligo.filt2$iCNV.malig.ident[which(names(oligo.filt2$orig.ident) %in%
                                     df.cell.types.thissamp$df.cell.types.cells)] = "non-malig"
oligo.filt2$iCNV.malig.ident[which(!(names(oligo.filt2$orig.ident) %in%
                                     df.cell.types.thissamp$df.cell.types.cells))] = "malignant"


# re-run initial processing, regressing
oligo.filt2 <- SCTransform(oligo.filt2, verbose = FALSE,
                     vars.to.regress = c('percent.mt', 'nCount_RNA'))
oligo.filt2 <- RunPCA(oligo.filt2)
oligo.filt2 <- FindNeighbors(oligo.filt2, dims = 1:30, verbose = FALSE)
oligo.filt2 <- FindClusters(oligo.filt2, verbose = FALSE, graph.name = "SCT_snn")
oligo.filt2 <- FindClusters(oligo.filt2, verbose = FALSE, graph.name = "SCT_snn", resolution = 1.2)
oligo.filt2 <- FindClusters(oligo.filt2, verbose = FALSE, graph.name = "SCT_snn", resolution = 0.3)
oligo.filt2 <- RunUMAP(oligo.filt2, dims = 1:30, reduction.key = "UMAP_", reduction.name = "umap.rna")

if("ATAC" %in% names(oligo.filt2@assays)){
  DefaultAssay(oligo.filt2) <- "ATAC"
  oligo.filt2 <- RunTFIDF(oligo.filt2)
  oligo.filt2 <- FindTopFeatures(oligo.filt2, min.cutoff = 'q0')
  oligo.filt2 <- RunSVD(oligo.filt2)
  oligo.filt2 <- FindNeighbors(object = oligo.filt2, reduction = 'lsi', dims = 2:50)
  oligo.filt2 <- FindClusters(object = oligo.filt2, verbose = FALSE, algorithm = 3, graph.name = "ATAC_snn")
  oligo.filt2 <- RunUMAP(oligo.filt2, reduction = 'lsi', dims = 2:50,
                   reduction.name = "umap.atac", reduction.key = "atacUMAP_")

  #compute WNN graph representing weighted combination of RNA and ATAC
  oligo.filt2 <- FindMultiModalNeighbors(oligo.filt2, reduction.list = list("pca", "lsi"),
                                   dims.list = list(1:50, 2:50))
  oligo.filt2 <- FindClusters(oligo.filt2, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  oligo.filt2 <- RunUMAP(oligo.filt2, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                   reduction.key = "wnnUMAP_")
}
saveRDS(oligo.filt2, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                             "/oligo_filt_final_", samples.newnames[pbs.arrayid], ".rds"))


# oligo.filt2 = readRDS(file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
#                                    "/oligo_filt_final_", samples.newnames[pbs.arrayid], ".rds"))

# #make malignant seurat object to save
# if(sum(oligo.filt2$iCNV.malig.ident=="malignant")>0){
#   DefaultAssay(oligo.filt2) = "RNA"
#   #get cells from clusters marked as malignant
#   clls = names(oligo.filt2$iCNV.malig.ident)[oligo.filt2$iCNV.malig.ident=="malignant"]
#   oligo.malig = subset(oligo.filt2, cells = clls)
# 
#   # normalize data with SCTransform()
#   oligo.malig <- SCTransform(oligo.malig, assay = 'RNA', new.assay.name = 'SCT',
#                              vars.to.regress = c('percent.mt', 'nCount_RNA'))
#   oligo.malig <- RunPCA(oligo.malig)
#   oligo.malig <- FindNeighbors(oligo.malig, dims = 1:30, verbose = FALSE)
#   oligo.malig <- FindClusters(oligo.malig, graph.name = "SCT_snn", verbose = FALSE)
#   oligo.malig = RunUMAP(oligo.malig, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',
#                         umap.method = "uwot", dims = 1:30)
# 
#   saveRDS(oligo.malig, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
#                         "/oligo_malig_", samples.newnames[pbs.arrayid], ".rds"))
# }

# if(sum(oligo.filt2$iCNV.malig.ident=="non-malig")>0){
#   DefaultAssay(oligo.filt2) = "RNA"
#   #get cells from clusters marked as malignant
#   clls = names(oligo.filt2$iCNV.malig.ident)[oligo.filt2$iCNV.malig.ident=="non-malig"]
#   oligo.malig = subset(oligo.filt2, cells = clls)
#   
#   # normalize data with SCTransform()
#   oligo.malig <- SCTransform(oligo.malig, assay = 'RNA', new.assay.name = 'SCT',
#                              vars.to.regress = c('percent.mt', 'nCount_RNA'))
#   oligo.malig <- RunPCA(oligo.malig)
#   oligo.malig <- FindNeighbors(oligo.malig, dims = 1:30, verbose = FALSE)
#   oligo.malig <- FindClusters(oligo.malig, graph.name = "SCT_snn", verbose = FALSE)
#   oligo.malig = RunUMAP(oligo.malig, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',
#                         umap.method = "uwot", dims = 1:30)
#   
#   saveRDS(oligo.malig, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
#                                      "/oligo_nonmalig_", samples.newnames[pbs.arrayid], ".rds"))
# }









##################
## make new reference and re-run inferCNV
##################



