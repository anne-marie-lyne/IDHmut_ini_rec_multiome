rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
script.dir = as.character(args[3])

#script to make non-malignant cell reference for inferCNV

data.dir = "/data/users/alyne/gliomult/svn/analyse/results"
samp.file = "/data/users/alyne/scrnaseq/oligo_10x_multiome/oligo_astro_samples_all.txt"
script.dir = "/data/users/alyne/scrnaseq/oligo_10x_multiome"

setwd(data.dir)

source(paste0(script.dir,"/make_cell_type_palette.R"))

library(Seurat)
library(Signac)

#number of cells for inferCNV reference
Nref = 8000

#reference used for SingleR
sr.ref = "imm"

#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
                         stringsAsFactors = F, fill = T)
samples = sample.info[,1]
samples.newnames = sample.info[,2]


#read in seurat objects for each sample
oligo.list = list()
for(i in 1:length(samples)){
  oligo.list[[i]] = readRDS(paste0(data.dir, "/analysis_", samples.newnames[i],"/oligo_filt_", samples.newnames[i], ".rds"))
  DefaultAssay(oligo.list[[i]]) = "RNA"
  oligo.list[[i]][['SCT']] <- NULL
  if("ATAC" %in% names(oligo.list[[i]]@assays)){
    oligo.list[[i]][['ATAC']] <- NULL
    oligo.list[[i]][["lsi"]] <- NULL
    oligo.list[[i]][["umap.atac"]] <- NULL
    oligo.list[[i]][["wnn.umap"]] <- NULL
  }
  oligo.list[[i]][["pca"]] <- NULL
  oligo.list[[i]][["umap.rna"]] <- NULL
}


#############################################
#### reference of Macrophage (myeloid) and Oligodendrocyte cells
#############################################

#combine cell types
if(sr.ref=="imm"){
  oligo.list = lapply(oligo.list, function(x) {
    ct.col = which(colnames(x@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
    x@meta.data[,ct.col] = combine_cell_types_imm(x@meta.data[,ct.col])
    return(x)
  })
}else{
  oligo.list = lapply(oligo.list, function(x) {
    ct.col = which(colnames(x@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
    x@meta.data[,ct.col] = combine_cell_types_Lin(x@meta.data[,ct.col])
    return(x)
  })
}


#count how many cells are macrophage or oligodendrocyte in each sample
types.nonmalig = c("Macrophage","Oligo")
Nnonmalig = sapply(oligo.list, function(x) {
  ct.col = which(colnames(x@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
  sum(x@meta.data[,ct.col] %in% types.nonmalig)
})
Ntotal = sum(Nnonmalig)

#number to sample from each sample
Nsample = sapply(Nnonmalig, function(x) {
  ceiling(x*Nref/Ntotal)
})
print(Nsample)
print(sum(Nsample))

#make list of filtered seurat objects with only oligo and mac cells
oligo.nonmalig.list = lapply(seq(length(oligo.list)), function(x) {
  ct.col = which(colnames(oligo.list[[x]]@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
  #find all cells which are non OPC or astrocyte
  clls = names(oligo.list[[x]]$orig.ident)[(oligo.list[[x]]@meta.data[,ct.col] %in% types.nonmalig)]

  if(Nsample[x]>0){
    #subsample the required number of cells
    clls.samp = clls[sample(length(clls),Nsample[x],replace = F)]
    #remove unnecessary assays
    oligo.samp = subset(oligo.list[[x]], cells = clls.samp)
    return(oligo.samp)
  }else{
    return(NA)
  }

})
oligo.nonmalig.list = oligo.nonmalig.list[!is.na(oligo.nonmalig.list)]

#get meta data for final Seurat object
cell.types = unlist(lapply(oligo.nonmalig.list, function(x) {
  ct.col = which(colnames(x@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
  ct = x@meta.data[,ct.col]
  names(ct) = paste(names(ct),x$orig.ident[1],sep="_")
  return(ct)
}))
orig.idents = unlist(lapply(oligo.nonmalig.list, function(x) {
  oi = x$orig.ident
  names(oi) = paste(names(oi),x$orig.ident[1],sep="_")
  return(oi)
}))
oligo.nonmalig.counts = lapply(oligo.nonmalig.list, function(x) {
  df1 = as.data.frame(GetAssayData(x, layer="counts"))
  colnames(df1) = paste(colnames(df1),x$orig.ident[1],sep="_")
  df2 = cbind(rownames(df1),df1)
  colnames(df2) = c("gene.names",colnames(df1))
  return(df2)
})

oligo.nonmalig.count = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene.names", all = TRUE), oligo.nonmalig.counts)
rownames(oligo.nonmalig.count) = oligo.nonmalig.count[,1]
oligo.nonmalig.count = oligo.nonmalig.count[,-1]

#make Seurat object
oligo.nonmalig = CreateSeuratObject(counts = oligo.nonmalig.count)
oligo.nonmalig@meta.data = cbind(oligo.nonmalig@meta.data,cell.types)
colnames(oligo.nonmalig@meta.data) = c(colnames(oligo.nonmalig@meta.data)[1:(length(colnames(oligo.nonmalig@meta.data))-1)],
                                       paste0('cell.type.sr.cl.',sr.ref))
oligo.nonmalig$orig.ident = orig.idents

#check if any cell types have fewer than two cells and remove if so
ct.col = which(colnames(oligo.nonmalig@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
Ntype = table(oligo.nonmalig@meta.data[,ct.col])
if(any(Ntype<2)){
  inds = which(Ntype<2)
  for(j in 1:length(inds)){
    inds2 = which(oligo.nonmalig@meta.data[,ct.col]==names(Ntype)[inds[j]])
    oligo.nonmalig = oligo.nonmalig[,-inds2]
  }
}

print(table(oligo.nonmalig@meta.data[,ct.col]))

#save reference
saveRDS(oligo.nonmalig, file = paste0(data.dir, "/oligo_nonmalig_merge_macrophage_oligo.rds"))



