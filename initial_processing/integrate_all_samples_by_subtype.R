#script to make figures summarising whole cohort

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
script.dir = as.character(args[3])
cancer = as.character(args[4])

library(Seurat)
library(Signac)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(SingleCellExperiment)
library(harmony)
library(readxl)
library(viridis)
library(UCell)
library(data.table)
library(future)
options(future.globals.maxSize = 120 * 1024 ^ 3)
pal = viridis(n=50)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


setwd(data.dir)
source(paste0(script.dir,"/make_cell_type_palette.R"))

sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
                         stringsAsFactors = F)
sample.info = sample.info[order(sample.info$V2),]
samples.newnames = sample.info[,2]
samples.newnames = sort(samples.newnames[grepl(cancer,samples.newnames,fixed=T)])


# remove the samples with no malignant cells
samples.newnames = samples.newnames[!(samples.newnames %in% c("oligo_P0003_rec2_S0007","astro_P0014_ini_S0026",
                             "oligo_P0022_rec1_S0043","oligo_P0030_ini_S0058"))]


st=1
en=length(samples.newnames)


#load all samples
#remove ATAC assay and as much stuff as possible to reduce size
oligo.filt.list = list()
for(i in st:en){
  print(i)
  oligo.filt.list[[i]] = readRDS(paste0(data.dir, "/analysis_", samples.newnames[i],
                                        "/oligo_filt_", samples.newnames[i], ".rds"))
  DefaultAssay(oligo.filt.list[[length(oligo.filt.list)]]) = "RNA"
  oligo.filt.list[[length(oligo.filt.list)]][["RNA"]] <- as(object = oligo.filt.list[[length(oligo.filt.list)]][["RNA"]], Class = "Assay5")
  oligo.filt.list[[length(oligo.filt.list)]][['SCT']] = NULL
  oligo.filt.list[[length(oligo.filt.list)]][["umap.rna"]] <- NULL
  oligo.filt.list[[length(oligo.filt.list)]][["pca"]] <- NULL
  if("ATAC" %in% names(oligo.filt.list[[length(oligo.filt.list)]]@assays)){
    oligo.filt.list[[length(oligo.filt.list)]][['ATAC']] = NULL
    oligo.filt.list[[length(oligo.filt.list)]][["lsi"]] <- NULL
    oligo.filt.list[[length(oligo.filt.list)]][["umap.atac"]] <- NULL
    oligo.filt.list[[length(oligo.filt.list)]][["wnn.umap"]] <- NULL
  }
  #rename cells
  oligo.filt.list[[length(oligo.filt.list)]] = RenameCells(oligo.filt.list[[length(oligo.filt.list)]],
                                     new.names = paste0(names(oligo.filt.list[[length(oligo.filt.list)]]$orig.ident),
                                                        "-",oligo.filt.list[[length(oligo.filt.list)]]$orig.ident[1]))
}

#add metadata to say whether sample is initial or recurrence and whether treated or untreated
oligo.filt.list = lapply(oligo.filt.list, function(x) {
  x$ini_rec = sapply(x$orig.ident, function(y) {
    unlist(strsplit(y,"_"))[3]
  })
  return(x)
})

comps <- as.data.frame(read_excel("sample_summary_astro_oligo_anon.xlsx", col_types = "text"))
oligo.filt.list = lapply(oligo.filt.list, function(x) {
  x$grade = as.numeric(comps$Grade[which(comps$Sample==x$orig.ident[1])])
  return(x)
})

oligo.filt.list = lapply(oligo.filt.list, function(x) {
  ind = as.numeric(comps$o_untreated_v_treated[which(comps$Sample==x$orig.ident[1])])
  if(is.na(ind)){
    ind = as.numeric(comps$a_untreated_v_treated[which(comps$Sample==x$orig.ident[1])])
  }
  x$treatment = c("untreated","treated")[ind]
  return(x)
})


#merge
oligo.ss = merge(oligo.filt.list[[1]], oligo.filt.list[2:length(oligo.filt.list)])

#remove doublet finder metadata as names are all slightly different
oligo.ss@meta.data = oligo.ss@meta.data[,!grepl("DF.classifications_", colnames(oligo.ss@meta.data), fixed=T)]
oligo.ss@meta.data = oligo.ss@meta.data[,!grepl("pANN_", colnames(oligo.ss@meta.data), fixed=T)]


#RPCA integration
oligo.ss <- SCTransform(oligo.ss, vars.to.regress = c('percent.mt', 'nCount_RNA'))
oligo.ss <- RunPCA(oligo.ss, npcs = 50, verbose = F)
oligo.ss <- IntegrateLayers(object = oligo.ss, method = RPCAIntegration, new.reduction = "integrated.rpca",
                            normalization.method = "SCT", verbose = F)
res = 0.3
oligo.ss <- FindNeighbors(oligo.ss, reduction = "integrated.rpca", dims = 1:50, graph.name = "rpca_nn")
oligo.ss <- FindClusters(oligo.ss, resolution = res, cluster.name = paste0("rpca_clusters",res), graph = "rpca_nn")
oligo.ss <- RunUMAP(oligo.ss, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")


#merge without integrating
DefaultAssay(oligo.ss) = "RNA"
layers_to_remove <- setdiff(Layers(oligo.ss[["RNA"]]), c("counts","data","scale.data"))
oligo.ss[["RNA"]]@layers <- oligo.ss[["RNA"]]@layers[!names(oligo.ss[["RNA"]]@layers) %in% layers_to_remove]
oligo.ss[["RNA"]] <- split(oligo.ss[["RNA"]], f = oligo.ss$orig.ident)

oligo.ss <- NormalizeData(oligo.ss)
oligo.ss <- FindVariableFeatures(oligo.ss)
oligo.ss <- ScaleData(oligo.ss)
oligo.ss <- RunPCA(oligo.ss)
res = 0.3
oligo.ss <- FindNeighbors(oligo.ss, dims = 1:50, reduction = "pca", graph.name = "unintegrated_nn")
oligo.ss <- FindClusters(oligo.ss, resolution = res, cluster.name = paste0("unintegrated_clusters",res), graph = "unintegrated_nn")
oligo.ss <- RunUMAP(oligo.ss, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")


#harmony integration
oligo.ss <- IntegrateLayers(object = oligo.ss, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.har", verbose = FALSE, assay = "RNA")
res = 0.3
oligo.ss <- FindNeighbors(oligo.ss, reduction = "integrated.har", dims = 1:50, graph.name = "har_nn")
oligo.ss <- FindClusters(oligo.ss, resolution = res, cluster.name = paste0("har_clusters",res), graph = "har_nn")
oligo.ss <- RunUMAP(oligo.ss, reduction = "integrated.har", dims = 1:50, reduction.name = "umap.har")


# Fmnn integration
library(SeuratWrappers)
oligo.ss <- IntegrateLayers(object = oligo.ss, method = FastMNNIntegration,
                            new.reduction = 'integrated.mnn', verbose = FALSE, batch = oligo.ss$orig.ident)
oligo.ss <- FindNeighbors(oligo.ss, reduction = "integrated.mnn", dims = 1:50, graph.name = "fmnn_nn")
res = 0.8
oligo.ss <- FindClusters(oligo.ss, resolution = res, cluster.name = paste0("fmnn_clusters",res), graph = "fmnn_nn")
oligo.ss <- RunUMAP(oligo.ss, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.fmnn")


if(cancer=="oligo"){
  #there are lots of clusters with small numbers of cells
  #move cells from clusters with fewer than 5 cells to cluster of nearest neighbour
  dim.str = "fmnn"
  res = 0.8
  md.col = which(colnames(oligo.ss@meta.data)==paste0(dim.str,"_clusters",res))
  f = oligo.ss@meta.data[,md.col]
  f = as.numeric(levels(f))[f]
  cl.count = table(f)
  cl.less5 = cl.count[cl.count<5]
  cl.greq5 = cl.count[cl.count>=5]
  nn_graph <- oligo.ss@graphs[[paste0(dim.str,"_nn")]]
  for(i in 1:length(cl.less5)){
    #find cells in this cluster
    cells = names(oligo.ss$orig.ident)[which(f==names(cl.less5)[i])]
    for(j in 1:length(cells)){
      neighbors <- names(which(nn_graph[cells[j], ] > 0))
      #find clusters of these cells
      cls = sapply(neighbors, function(x) {f[which(names(oligo.ss$orig.ident)==x)]})
      cl.to.assign = names(sort(table(cls), decreasing=T))[names(sort(table(cls), decreasing=T)) %in% names(cl.greq5)][1]
      oligo.ss@meta.data[which(names(oligo.ss$orig.ident)==cells[j]),md.col] = as.numeric(cl.to.assign)
    }
  }
  oligo.ss@meta.data[,md.col] = droplevels(oligo.ss@meta.data[,md.col])
  level.to.add = max(as.numeric(levels(oligo.ss@meta.data[,md.col])))+1
  levels(oligo.ss@meta.data[,md.col]) = c(levels(oligo.ss@meta.data[,md.col]),level.to.add)
  
  #split cluster 16 in 2 as these seem to be endothelial and T-cells mixed together
  oligo.ss@meta.data[which(oligo.ss@meta.data[,md.col]==16 & Embeddings(oligo.ss,"umap.fmnn")[,1]>(-4)),md.col] = level.to.add
}


oligo.ss$IDH.mut[which(oligo.ss$IDH.mut=="non_mutant" | oligo.ss$IDH.mut=="not_expressed")] = NA
oligo.ss$IDH.mut = factor(oligo.ss$IDH.mut)


#make atlas label ignoring adult/fetal to make list shorter
oligo.ss$atlas.label2 = sapply(oligo.ss$atlas.label, function(x) {
  unlist(strsplit(x,"_"))[1]
})


saveRDS(oligo.ss, paste0("oligo.ss_integrated_RPCA_all_",cancer,".rds"))




oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_all_",cancer,".rds"))


DefaultAssay(oligo.ss) ="RNA"
oligo.ss = JoinLayers(oligo.ss)


#plot umaps of all cells showing treatment status, ini_rec, cell type, grade etc


#also plot umap showing all cells with SNVs detected
#make data frame so can plot cells using ggplot with SNVs on top
dim.red = "umap.fmnn"
dim.str = "fmnn"
if(cancer=="astro"){
  res = 0.3
}else{
  res = 0.8
}
md.col = which(colnames(oligo.ss@meta.data)==paste0(dim.str,"_clusters",res))
if(cancer=="oligo"){
  df = data.frame(umap_1 = Embeddings(oligo.ss,dim.red)[,1], umap_2 = Embeddings(oligo.ss,dim.red)[,2], SNV = oligo.ss$IDH.mut,
                  cluster = oligo.ss@meta.data[,md.col], p1.ATAC = oligo.ss$p1.av.ATAC, q19.ATAC = oligo.ss$q19.av.ATAC)
}else{
  df = data.frame(umap_1 = Embeddings(oligo.ss,dim.red)[,1], umap_2 = Embeddings(oligo.ss,dim.red)[,2], SNV = oligo.ss$IDH.mut,
                  cluster = oligo.ss@meta.data[,md.col], q7.ATAC = oligo.ss$q7.av.ATAC)
}
levs = unique(df$SNV)[!is.na(unique(df$SNV))]
cols = gg_color_hue(length(levs))
df$color = sapply(df$SNV, function(x) {
  if(is.na(x)){return("gray70")}else{cols[which(levs==x)]}
})
#sort to put all NAs first, i.e., to plot all the other things last
df = df[order(df$SNV, na.last = F),]


pt.sz = 0.001
p1 = DimPlot(oligo.ss, group.by = "orig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p2 = DimPlot(oligo.ss, group.by = "cell.type.sr.cl.imm", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
if(cancer=="astro"){
  p3 = FeaturePlot(oligo.ss, features = "q7.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p4 = FeaturePlot(oligo.ss, features = "q7.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p5 = FeaturePlot(oligo.ss, features = "altchr.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
}else{
  p3 = FeaturePlot(oligo.ss, features = "p1.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p4 = FeaturePlot(oligo.ss, features = "p1.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p5 = FeaturePlot(oligo.ss, features = "q19.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
}
p6 = DimPlot(oligo.ss, group.by = "nCNA", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p7 = DimPlot(oligo.ss, group.by = "treatment", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p8 = DimPlot(oligo.ss, group.by = "ini_rec", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p9 = DimPlot(oligo.ss, group.by = "grade", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p10 = DimPlot(oligo.ss, group.by = "cell.type.sr.imm", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p11 = ggplot(df, aes(x = umap_1, y = umap_2, color = SNV)) + geom_point()
p12 = DimPlot(oligo.ss, group.by = paste0(dim.str,"_clusters",res), reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p13 = DimPlot(oligo.ss, group.by = "iCNV.malig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p14 = DimPlot(oligo.ss, group.by = "atlas.label2", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
if(cancer=="oligo"){
  p15 = ggplot(df, aes(x=cluster, y=p1.ATAC, fill=cluster)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0))
  p16 = ggplot(df, aes(x=cluster, y=q19.ATAC, fill=cluster)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0))
}else{
  p15 = ggplot(df, aes(x=cluster, y=q7.ATAC, fill=cluster)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0))
  p16 = plot_spacer()
}
print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p15 + p16 + plot_layout(ncol = 4))
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_",dim.str,".pdf"), width=24, height=19)



#find markers of each cluster versus all others
dim.red = "umap.fmnn"
dim.str = "fmnn"
if(cancer=="astro"){
  res = 0.3
}else{
  res = 0.8
}
# oligo.ss = PrepSCTFindMarkers(oligo.ss)
oligo.ss <- JoinLayers(oligo.ss)
all.markers <- FindAllMarkers(object = oligo.ss, assay = "RNA", slot = "data", group.by = paste0("fmnn_clusters",res))
fwrite(all.markers, file = paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_",dim.str,".tsv.gz"), sep="\t")



#make plots for each type of integration with clusters
pt.sz = 0.002
count = 1
plot.list = list()
for(dim.red in c("umap.unintegrated","umap.fmnn","umap.rpca","umap.har")){
  if(dim.red == "umap.fmnn"){
    dim.str = "fmnn"
    reduc = "integrated.mnn"
  }else if(dim.red == "umap.rpca"){
    dim.str = "rpca"
    reduc = "integrated.rpca"
  }else if(dim.red == "umap.unintegrated"){
    dim.str = "unintegrated"
    reduc = "integrated.dr"
  }else if(dim.red == "umap.har"){
    dim.str = "har"
    reduc = "integrated.har"
  }
  plot.list[[count]] = DimPlot(oligo.ss, group.by = paste0(dim.str,"_clusters",res), reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
  plot.list[[count]] = DimPlot(oligo.ss, group.by = "cell.type.sr.imm", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
  plot.list[[count]] = DimPlot(oligo.ss, group.by = "iCNV.malig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
}
wrap_plots(plot.list, ncol = 3)
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_compare_dim_red.pdf"), width=19, height=16)



#plot markers to confirm doublets
markers = c("KCNN3","FERMT1","PCDH15","TOP2A","PTPRC","CD3G","PLP1","LAMA2","SKAP1","RBFOX3","SLC17A7","RBPMS")
FeaturePlot(oligo.ss, features = markers, ncol = 3, reduction = "umap.fmnn", pt.size = pt.sz, raster = F)
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_markers_fmnn.pdf"), width=19, height=16)


#UP TO HERE


#adjustments to make to annotations and doublets
oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_all_",cancer,".rds"))

if(cancer=="astro"){
  res = 0.3
}else{
  res = 0.8
}
cl.use = "fmnn_clusters"
dim.red = "umap.fmnn"
md.col = which(colnames(oligo.ss@meta.data)==paste0(cl.use,res))

df.cell.types = data.frame(cell = names(oligo.ss$orig.ident), cluster = oligo.ss@meta.data[,md.col],
                           cluster_rpca = oligo.ss$rpca_clusters0.3, umap_1 = Embeddings(oligo.ss,dim.red)[,1],
                           umap_2 = Embeddings(oligo.ss,dim.red)[,2])
if(cancer=="astro"){
  oligo.ss$cell.type.sr.cl.imm.redef = oligo.ss$cell.type.sr.cl.imm
  #remove cells which seem to be doublets, i.e. expressing markers of two cell types and placed inbetween two cell types
  #cl 11 = malignant macrophage
  #cl 14 = malignant oligodendrocyte 
  #cl at the bottom which seems to express markers of everything
  cells.remove = names(oligo.ss$orig.ident)[which(oligo.ss@meta.data[,md.col]==11 | oligo.ss@meta.data[,md.col]==14)]
  cells.remove = c(cells.remove, rownames(df.cell.types)[which(abs(df.cell.types$umap_1)<2.5 & df.cell.types$umap_2<(-11))])
  cells.remove = unique(cells.remove)
  
  oligo.ss$cell.type.sr.cl.imm.redef[which(names(oligo.ss$orig.ident) %in% cells.remove)] = "DOUBLET"
  
  #redefine TME cell types based on fmn clusters
  cl.to.relabel = c(10,13,15,16,4,12,7)
  names.relabel = c("Glutamatergic","GABAergic","Astro","T-cells","Macrophage","Macrophage","Oligo")
  endo.cl = 11
  df.cell.types = df.cell.types[which((df.cell.types$cluster %in% cl.to.relabel) | df.cell.types$cluster_rpca==endo.cl),]
  df.cell.types$cell.type = rep("",nrow(df.cell.types))
  for(j in 1:length(cl.to.relabel)){
    df.cell.types$cell.type[which(df.cell.types$cluster==cl.to.relabel[j])] = names.relabel[j]
    oligo.ss$cell.type.sr.cl.imm.redef[which(oligo.ss@meta.data[,md.col]==cl.to.relabel[j])] = names.relabel[j]
  }
  
  #add small endothelial cluster from RPCA
  df.cell.types$cell.type[which(df.cell.types$cluster_rpca==endo.cl)] = "Endo"
  oligo.ss$cell.type.sr.cl.imm.redef[which(oligo.ss$rpca_clusters0.3==endo.cl)] = "Endo"

  #ensure this matches the malignant identity
  oligo.ss$iCNV.malig.ident[which(oligo.ss@meta.data[,md.col] %in% cl.to.relabel)] = "non−malig"
  oligo.ss$iCNV.malig.ident[which(!(oligo.ss@meta.data[,md.col] %in% cl.to.relabel))] = "malignant"
  oligo.ss$iCNV.malig.ident[which(oligo.ss$rpca_clusters0.3==endo.cl)] = "non−malig"
  
  oligo.ss$iCNV.malig.ident[which(oligo.ss$cell.type.sr.cl.imm.redef=="DOUBLET")] = "DOUBLET"
  
}else if(cancer=="oligo"){
  oligo.ss$cell.type.sr.cl.imm.redef = oligo.ss$cell.type.sr.cl.imm
  #remove cells which seem to be doublets, i.e. expressing markers of two cell types and placed inbetween two cell types
  #cl 13 = malignant macrophage
  #cl 12 = oligodendrocyte malignant
  #cl 18 = macrophage oligodendrocyte
  cells.remove = names(oligo.ss$orig.ident)[which(oligo.ss@meta.data[,md.col]==13 | oligo.ss@meta.data[,md.col]==12 |
                                                    oligo.ss@meta.data[,md.col]==18)]
  cells.remove = unique(cells.remove)
  
  oligo.ss$cell.type.sr.cl.imm.redef[which(names(oligo.ss$orig.ident) %in% cells.remove)] = "DOUBLET"
  
  #redefine TME cell types based on fmn clusters
  cl.to.relabel = c(16,14,7,19,5,17,11,15)
  names.relabel = c("T-cells","Macrophage","Macrophage","Endo","Oligo","Astro","Glutamatergic","GABAergic")
  df.cell.types = df.cell.types[which(df.cell.types$cluster %in% cl.to.relabel),]
  df.cell.types$cell.type = rep("",nrow(df.cell.types))
  for(j in 1:length(cl.to.relabel)){
    df.cell.types$cell.type[which(df.cell.types$cluster==cl.to.relabel[j])] = names.relabel[j]
    oligo.ss$cell.type.sr.cl.imm.redef[which(oligo.ss@meta.data[,md.col]==cl.to.relabel[j])] = names.relabel[j]
  }
  
  #ensure this matches the malignant identity
  oligo.ss$iCNV.malig.ident[which(oligo.ss@meta.data[,md.col] %in% cl.to.relabel)] = "non−malig"
  oligo.ss$iCNV.malig.ident[which(!(oligo.ss@meta.data[,md.col] %in% cl.to.relabel))] = "malignant"
  
  oligo.ss$iCNV.malig.ident[which(oligo.ss$cell.type.sr.cl.imm.redef=="DOUBLET")] = "DOUBLET"
}


saveRDS(oligo.ss, paste0("oligo.ss_integrated_RPCA_all_",cancer,".rds"))


#write file of probable doublets
write.table(names(oligo.ss$orig.ident)[which((names(oligo.ss$orig.ident) %in% cells.remove))],
            file = paste0("probable_doublets_",cancer,".txt"), col.names=F)
write.table(df.cell.types,
            file = paste0("cell_types_after_merge_",cancer,".txt"), col.names=F)


#replot figure using all integration methods
pt.sz = 0.002
count = 1
plot.list = list()
for(dim.red in c("umap.fmnn","umap.rpca","umap.har")){
  res = 0.3
  if(dim.red == "umap.fmnn"){
    dim.str = "fmnn"
    reduc = "integrated.mnn"
    if(cancer=="oligo"){
      res = 0.8
    }
  }else if(dim.red == "umap.rpca"){
    dim.str = "rpca"
    reduc = "integrated.rpca"
  }else if(dim.red == "umap.unintegrated"){
    dim.str = "unintegrated"
    reduc = "integrated.dr"
  }else if(dim.red == "umap.har"){
    dim.str = "har"
    reduc = "integrated.har"
  }
  plot.list[[count]] = DimPlot(oligo.ss, group.by = paste0(dim.str,"_clusters",res), reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
  plot.list[[count]] = DimPlot(oligo.ss, group.by = "cell.type.sr.cl.imm.redef", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
  plot.list[[count]] = DimPlot(oligo.ss, group.by = "iCNV.malig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
  count = count + 1
}
wrap_plots(plot.list, ncol = 3)
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_compare_dim_red_redef.pdf"), width=19, height=12)


DimPlot(oligo.ss, group.by = "cell.type.malig", reduction = dim.red, split.by = "orig.ident", raster = F, ncol=6)
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_by_sample.pdf"), width=19, height=12)
#
dim.red = "umap.fmnn"
oligo.ss$cell.type.malig = oligo.ss$cell.type.sr.cl.imm.redef
oligo.ss$cell.type.malig[which(oligo.ss$iCNV.malig.ident=="malignant")] = "malignant"
plot.list = list()
for(i in 1:length(unique(oligo.ss$orig.ident))){
  oligo.subset = subset(oligo.ss, subset = orig.ident==unique(oligo.ss$orig.ident)[i])
  plot.list[[i]] = DimPlot(oligo.subset, group.by = "cell.type.malig", reduction = dim.red) + 
    ggtitle(unique(oligo.ss$orig.ident)[i])
}
wrap_plots(plot.list)
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_by_sample.pdf"), width=19, height=12)


#make full plot with chosen dim reduc
oligo.ss = subset(oligo.ss, cells = names(oligo.ss$orig.ident)[which(!(names(oligo.ss$orig.ident) %in% cells.remove))])
saveRDS(oligo.ss, paste0("oligo.ss_integrated_RPCA_all_",cancer,"_doub_rem.rds"))


res = 0.8
cl.use = "fmnn_clusters"
dim.red = "umap.fmnn"
dim.str = "fmnn"
md.col = which(colnames(oligo.ss@meta.data)==paste0(cl.use,res))

#make data frame so can plot SNVs last
df = data.frame(umap_1 = Embeddings(oligo.ss,dim.red)[,1],
                umap_2 = Embeddings(oligo.ss,dim.red)[,2], SNV = oligo.ss$IDH.mut)
levs = unique(df$SNV)[!is.na(unique(df$SNV))]
cols = gg_color_hue(length(levs))
df$color = sapply(df$SNV, function(x) {
  if(is.na(x)){return("gray70")}else{cols[which(levs==x)]}
})
#sort to put all NAs first
df = df[order(df$SNV, na.last = F),]

pt.sz = 0.002
p1 = DimPlot(oligo.ss, group.by = "orig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p2 = DimPlot(oligo.ss, group.by = "cell.type.sr.cl.imm.redef", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
if(cancer=="astro"){
  p3 = FeaturePlot(oligo.ss, features = "q7.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p4 = FeaturePlot(oligo.ss, features = "q7.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p5 = FeaturePlot(oligo.ss, features = "altchr.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
}else{
  p3 = FeaturePlot(oligo.ss, features = "p1.av", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p4 = FeaturePlot(oligo.ss, features = "p1.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
  p5 = FeaturePlot(oligo.ss, features = "q19.av.ATAC", reduction = dim.red, pt.size = pt.sz, raster=F) +
    scale_colour_gradientn(colours = pal, na.value = "grey90")
}
p6 = DimPlot(oligo.ss, group.by = "nCNA", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p7 = DimPlot(oligo.ss, group.by = "treatment", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p8 = DimPlot(oligo.ss, group.by = "ini_rec", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p9 = DimPlot(oligo.ss, group.by = "grade", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p10 = DimPlot(oligo.ss, group.by = "cell.type.sr.imm", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p11 = ggplot(df, aes(x = umap_1, y = umap_2, color = SNV)) + geom_point()
p12 = DimPlot(oligo.ss, group.by = paste0(cl.use,res), reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
p13 = DimPlot(oligo.ss, group.by = "iCNV.malig.ident", reduction = dim.red, pt.size = pt.sz, shuffle=T, raster=F)
print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + plot_layout(ncol = 4))
ggsave(paste0(data.dir,"/figures/umap_ss_all_cells_",cancer,"_",dim.str,"_doubrem.pdf"), width=24, height=16)






# #make figure with expression of known marker genes
# #and boxplot of ATAC scores for 1p and 19q
# mark.genes = c("AQP4","APOE","SLC2A1",#astrocytes
#                "CSPG4","OLIG2","PDGFRA",#OPC
#                "MBP","MAG","MOG",#oligodendrocytes
#                "RBFOX3","MAP2","GAP43",#neurons
#                "SLC17A7","SLC17A6","GRIN1",#glutamatergic
#                "SLC6A1","GAD1","GAD2",#gabaergic
#                "PTPRC","CD3E","CD3D",#T-cells
#                "CCL3", "ITGAM", "TMEM119",#microglia
#                "CD40LG","HLA-DRB1","CD19",#B-cells
#                "FLT1", "CLDN5", "ABCB1",#endothelial
#                "OGN","LUM","DCN"#VLMC
#                )
