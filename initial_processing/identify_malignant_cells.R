
rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
pbs_arrayid = as.numeric(args[3])
ref = as.character(args[4])
hmm = as.logical(args[5])
script.dir = as.character(args[6])

#script to plot lots of different information about cells in each sample to identify malignant cells
#SNVs detected, inferCNV score etc


library(data.table)
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(R.utils)
library(RColorBrewer)
library(viridis)

source(paste0(script.dir,"/make_cell_type_palette.R"))
source(paste0(script.dir,"/define_malignant_function.R"))

#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
                         stringsAsFactors = F, fill = T)
samples = sample.info[,1]
samples.newnames = sample.info[,2]
samp = samples.newnames[pbs_arrayid]
print(samp)

#load Seurat object
oligo.filt = readRDS(paste0(data.dir,"/analysis_",samp,"/oligo_filt_",samp,".rds"))


#load inferCNV reference
oligo.nonmalig = readRDS(paste0(data.dir,"/oligo_nonmalig_merge_",ref,".rds"))
cl.res = 0.8 #clustering resolution
md.col = which(colnames(oligo.filt@meta.data)==paste0("SCT_snn_res.",cl.res))
print(md.col)

#amount above or below 1 infercnv score needs to be to be defined malignant
icnv.diff = 0.025

pal = viridis(n=50)

#reference used for SingleR
sr.ref = "imm"


#plot umpa with clusters
p1 = DimPlot(oligo.filt, group.by = paste0('SCT_snn_res.',cl.res), reduction = 'umap.rna', pt.size = 0.005) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
ct.col = which(colnames(oligo.filt@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
ct.col.nm = which(colnames(oligo.nonmalig@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
oligo.filt@meta.data[,ct.col] = combine_cell_types_imm(oligo.filt@meta.data[,ct.col])


#plot umap with SingleR cell type
p2 = DimPlot(oligo.filt, group.by = paste0('cell.type.sr.cl.',sr.ref), reduction = 'umap.rna', pt.size = 0.005) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")

############################
#### rough count of number CNA from infercnv output
############################

#input reference inferCNV
ref.pth = paste0(data.dir,"/analysis_",samp,"/inferCNV_",samp,"_hmm_",as.character(hmm),"_",ref)
iCNV.ref = fread(paste0(ref.pth,"/infercnv.references.txt.gz"),data.table = F)
rownames(iCNV.ref) = iCNV.ref[,1]
iCNV.ref = iCNV.ref[,2:ncol(iCNV.ref)]
# #take only cells from this sample
# iCNV.ref = iCNV.ref[,grepl(samp,colnames(iCNV.ref),fixed=T)]

#input obs inferCNV
iCNV.obs = fread(paste0(ref.pth,"/infercnv.observations.txt.gz"),data.table = F)
rownames(iCNV.obs) = iCNV.obs[,1]
iCNV.obs = iCNV.obs[,2:ncol(iCNV.obs)]

#input gene order file
gene.order = fread(paste0(data.dir,"/gene_order_file.txt"),data.table = F)

#compute average for each gene across all test cells
gene.av = rowMeans(iCNV.obs)

#make plot of CNA v genome position
df = data.frame(gene.pos = seq(length(gene.av)), gene.av = gene.av)
ggplot(df, aes(gene.pos, gene.av)) + geom_point(size = 2)
ggsave(paste0(data.dir, "/figures_infercnv/CNAs_average_",samp,".pdf"),
       width = 10, height = 4)

#find genes with average scores over 1.01 or under 0.99 
#(quite weak to pick up subclonal)
inds = which(gene.av>1.01 | gene.av<0.99)

#find separate regions
region.list = list()
j=1
region.list[[j]] = inds[1]
gene1 = c(rownames(iCNV.obs)[inds[1]])
for(i in 2:length(inds)){
  if(abs(inds[i]-inds[i-1])==1){
    region.list[[j]] = c(region.list[[j]],inds[i])
  }else{
    j=j+1
    region.list[[j]] = inds[i]
    gene1 = c(gene1, rownames(iCNV.obs)[inds[i]])
  }
}
#keep only CNA regions at least 150 genes long
min.length = 150
gene1 = gene1[sapply(region.list, function(x) {length(x)>=min.length})]
region.list = region.list[sapply(region.list, function(x) {length(x)>=min.length})]
print(sapply(region.list, length))
#find which chromosome first gene is on
chroms = sapply(gene1, function(x) {gene.order[which(gene.order[,1]==x),2]})

#if sample has some CNA
if(length(region.list)>0){
  #compute average score for each region for each cell
  ref.scores = lapply(region.list, function(x) {
    colMeans(iCNV.ref[x,])
  })
  obs.scores = lapply(region.list, function(x) {
    colMeans(iCNV.obs[x,])
  })
  
  #make boxplot of scores by cluster/cell type for each CNA region
  plot.list = lapply(seq(length(region.list)), function(y) {
    #get scores in test clusters
    scores.obs.cl = unlist(lapply(unique(oligo.filt@meta.data[,md.col]), function(x) {
      #find all cells from that cluster
      clls = names(oligo.filt$orig.ident)[which(oligo.filt@meta.data[,md.col]==x)]
      obs.scores[[y]][sapply(names(obs.scores[[y]]), 
                             function(z) {unlist(strsplit(z,"_",fixed=T))[1]}) %in% clls]
    }))
    #get scores in reference cells
    scores.ref.type = unlist(lapply(unique(oligo.nonmalig@meta.data[,ct.col.nm]), function(x) {
      #find all cells from that cell type
      clls = names(oligo.nonmalig@meta.data[,ct.col.nm])[which(oligo.nonmalig@meta.data[,ct.col.nm]==x)]
      ref.scores[[y]][sapply(names(ref.scores[[y]]), function(z) {z %in% clls})]
    }))
    scores = c(scores.obs.cl, scores.ref.type)
    clusters.obs.cl = unlist(lapply(unique(oligo.filt@meta.data[,md.col]), function(x) {
      clls = names(oligo.filt$orig.ident)[which(oligo.filt@meta.data[,md.col]==x)]
      rep(as.numeric(paste(x)), sum(sapply(names(obs.scores[[y]]), 
                                           function(z) {unlist(strsplit(z,"_",fixed=T))[1]}) %in% clls))
    }))
    clusters.ref.type = unlist(lapply(unique(oligo.nonmalig@meta.data[,ct.col.nm]), function(x) {
      clls = names(oligo.nonmalig@meta.data[,ct.col.nm])[which(oligo.nonmalig@meta.data[,ct.col.nm]==x)]
      rep(x, sum(sapply(names(ref.scores[[y]]), function(z) {z %in% clls})))
    }))
    clusters = c(clusters.obs.cl,clusters.ref.type)
    df = data.frame(score = scores, cluster = clusters)
    ggplot(df, aes(x=cluster, y=score)) + geom_boxplot() + ggtitle(chroms[y]) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  })
  p1.list = wrap_plots(plot.list)
  print(p1.list)
  ggsave(paste0(data.dir, "/figures_infercnv/boxplots_iCNV_region_scores_",samp,".png"), 
         width = 15, height = 10)
  
  
  #find quantiles for each cluster
  region.scores = lapply(seq(length(region.list)), function(y) {
    lapply(unique(oligo.filt@meta.data[,md.col]), function(x) {
      #find all cells from that cluster
      clls = names(oligo.filt$orig.ident)[which(oligo.filt@meta.data[,md.col]==x)]
      #compute median and quantile values for all cells in that cluster
      quantile(c(ref.scores[[y]][sapply(names(ref.scores[[y]]), 
                                        function(z) {unlist(strsplit(z,"_",fixed=T))[1]}) %in% clls],
                 obs.scores[[y]][sapply(names(obs.scores[[y]]), 
                                        function(z) {unlist(strsplit(z,"_",fixed=T))[1]}) %in% clls]),
               probs = c(0.25,0.5,0.75))
    })
  })
  
  #for each cluster, count how many CNAs it has
  #for gains: lower quartile has to be >1 and median >1.02; for losses: upper quartile has to be <1 and median <0.98
  nCNA.clus = sapply(seq(length(region.scores[[1]])), function(x) {
    sum(sapply(region.scores, function(y) {
      (y[[x]][1]>1 && y[[x]][2]>1.02) || (y[[x]][3]<1 && y[[x]][2]<0.98)
    }))
  })
  
  #compute value for each cell
  nCNA = sapply(oligo.filt@meta.data[,md.col], function(x) {
    nCNA.clus[which(unique(oligo.filt@meta.data[,md.col])==x)]
  })
}else{
  nCNA = rep(0,length(oligo.filt$orig.ident))
}
names(nCNA) = names(oligo.filt$orig.ident)
oligo.filt$nCNA = nCNA

#plot umap with rough count of CNA for each test cell and for cells which were in reference
p3 = DimPlot(oligo.filt, group.by = 'nCNA', reduction = 'umap.rna', pt.size = 0.005) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")



###########################################
#### plot average of 1p and 7q genes from RNA infercnv 
###########################################

#locations of p/q arms
pq.loc = fread(paste0(data.dir, "/cytoBand.txt.gz"), data.table = F)
p1.st = 0 #start of 1p
p1.en = pq.loc[which(pq.loc[,1]=="chr1" & grepl("^q",pq.loc[,4]))[1]-1,3] #end of 1p
#get all genes in this region on 1p
p1.genes = gene.order[gene.order[,2]=="chr1" & gene.order[,4]<p1.en,1]
#average inferCNV for 1p genes in reference cells
p1.ref.av = colMeans(iCNV.ref[rownames(iCNV.ref) %in% p1.genes,])

#add 1p averages to seurat object for each cell
cl.uni = unique(oligo.filt@meta.data[,md.col])

#average for observation cells
p1.obs.av = colMeans(iCNV.obs[rownames(iCNV.obs) %in% p1.genes,])
p1.av = sapply(paste(colnames(oligo.filt@assays$SCT), oligo.filt$orig.ident[1], sep = "_"),
               function(y) {
                 if(y %in% names(c(p1.obs.av,p1.ref.av))){
                   c(p1.obs.av,p1.ref.av)[which(names(c(p1.obs.av,p1.ref.av))==y)]
                 }else{
                   NA
                 }
               })
names(p1.av) = colnames(oligo.filt@assays$SCT)
oligo.filt$p1.av = p1.av


#same for 7q amplification
q7.rows = which(pq.loc[,1]=="chr7" & grepl("^q",pq.loc[,4]))
q7.st = pq.loc[q7.rows[1],2] #start of 7q
q7.en = pq.loc[q7.rows[length(q7.rows)],3] #end of 7q
#get all genes in this region on 7q
q7.genes = gene.order[gene.order[,2]=="chr7" & gene.order[,3]>q7.st & gene.order[,4]<q7.en,1]
#average inferCNV for 7q genes in reference cells
#try running using only last 200 genes to get cleaner signal
logi = rownames(iCNV.ref) %in% q7.genes
#find number of 7q genes
n7q = sum(logi)
inds7q = which(logi)
#leave only last 220 TRUE, set all others to FALSE
logi[inds7q[1:(n7q-220)]] = FALSE
q7.ref.av = colMeans(iCNV.ref[logi,])

#add q7 averages to seurat object for each cell
#average for observation cells
q7.obs.av = colMeans(iCNV.obs[logi,])
q7.av = sapply(paste(colnames(oligo.filt@assays$SCT),
                     oligo.filt$orig.ident[1], sep = "_"),
               function(y) {
                 if(y %in% names(c(q7.obs.av,q7.ref.av))){
                   c(q7.obs.av,q7.ref.av)[which(names(c(q7.obs.av,q7.ref.av))==y)]
                 }else{
                   NA
                 }
               })
names(q7.av) = colnames(oligo.filt@assays$SCT)
oligo.filt$q7.av = q7.av

#find doublet score column to add to data frame for figure later
md.doub = grepl("pANN_",colnames(oligo.filt@meta.data),fixed=T)


#plot umap and boxplot by cluster of 1p/7q mean for oligo/astro respectively
if(grepl("oligo",samp)){
  p4 = FeaturePlot(oligo.filt, features = 'p1.av', reduction = 'umap.rna', pt.size = 0.005) +
    scale_colour_gradientn(colours = pal, na.value = "grey90") +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  
  #boxplot showing 1p per cluster
  df = data.frame(cluster = oligo.filt@meta.data[,md.col], p1 = oligo.filt$p1.av,
                  cell.type = oligo.filt@meta.data[,ct.col], doub.score = oligo.filt@meta.data[,md.doub])
  p5 = ggplot(df, aes(x=cluster, y=p1, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
    geom_hline(yintercept=(1-icnv.diff), color = "blue", linewidth=1.5)
}else{
  p4 = FeaturePlot(oligo.filt, features = 'q7.av', reduction = 'umap.rna', pt.size = 0.005) +
    scale_colour_gradientn(colours = pal, na.value = "grey90") +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  df = data.frame(cluster = oligo.filt@meta.data[,md.col], q7 = oligo.filt$q7.av,
                  cell.type = oligo.filt@meta.data[,ct.col], doub.score = oligo.filt@meta.data[,md.doub])
  p5 = ggplot(df, aes(x=cluster, y=q7, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
    geom_hline(yintercept=(1+icnv.diff), color = "blue", linewidth=1.5) +
    coord_cartesian(ylim = c(0.95,1.13))
}



###########################################
#### plot average of 1p and 7q genes from ATAC infercnv if ATAC assay exists
###########################################

ref.pth = paste0(data.dir,"/analysis_",samp,"/infercnv_ATAC_",samp)

if(file.exists(paste0(ref.pth,"/infercnv.references.txt.gz"))){
  #input reference inferCNV
  iCNV.ref.ATAC = fread(paste0(ref.pth,"/infercnv.references.txt.gz"),data.table = F)
  rownames(iCNV.ref.ATAC) = iCNV.ref.ATAC[,1]
  iCNV.ref.ATAC = iCNV.ref.ATAC[,2:ncol(iCNV.ref.ATAC)]
  #input obs inferCNV
  iCNV.obs.ATAC = fread(paste0(ref.pth,"/infercnv.observations.txt.gz"),data.table = F)
  rownames(iCNV.obs.ATAC) = iCNV.obs.ATAC[,1]
  iCNV.obs.ATAC = iCNV.obs.ATAC[,2:ncol(iCNV.obs.ATAC)]
  #input gene order file
  gene.order.ATAC = fread(paste0(ref.pth,"/",samp,"_ATAC_cnv_ref.txt.gz"),data.table = F)
  #locations of p/q arms
  pq.loc = fread(paste0(data.dir, "/cytoBand.txt.gz"), data.table = F)
  p1.st = 0 #start of 1p
  p1.en = pq.loc[which(pq.loc[,1]=="chr1" & grepl("^q",pq.loc[,4]))[1]-1,3] #end of 1p
  #get all genes in this region on 1p
  p1.genes = gene.order.ATAC[gene.order.ATAC[,2]=="1" & gene.order.ATAC[,4]<p1.en,1]
  #average inferCNV for 1p genes in reference cells
  p1.ref.av = colMeans(iCNV.ref.ATAC[rownames(iCNV.ref.ATAC) %in% p1.genes,])
  #add 1p averages to seurat object for each cell
  cl.uni = unique(oligo.filt@meta.data[,md.col])
  #average for observation cells
  p1.obs.av = colMeans(iCNV.obs.ATAC[rownames(iCNV.obs.ATAC) %in% p1.genes,])
  p1.av = sapply(colnames(oligo.filt@assays$SCT),
                 function(y) {
                   if(y %in% names(c(p1.obs.av,p1.ref.av))){
                     c(p1.obs.av,p1.ref.av)[which(names(c(p1.obs.av,p1.ref.av))==y)]
                   }else{
                     NA
                   }
                 })
  names(p1.av) = colnames(oligo.filt@assays$SCT)
  oligo.filt$p1.av.ATAC = p1.av
  #same for 7q amplification
  q7.rows = which(pq.loc[,1]=="chr7" & grepl("^q",pq.loc[,4]))
  q7.st = pq.loc[q7.rows[1],2] #start of 7q
  q7.en = pq.loc[q7.rows[length(q7.rows)],3] #end of 7q
  #get all genes in this region on 7q
  q7.genes = gene.order.ATAC[gene.order.ATAC[,2]=="7" & gene.order.ATAC[,3]>q7.st & gene.order.ATAC[,4]<q7.en,1]
  #average inferCNV for 7q genes in reference cells
  logi = rownames(iCNV.ref.ATAC) %in% q7.genes
  #find number of 7q genes
  n7q = sum(logi)
  inds7q = which(logi)
  q7.ref.av = colMeans(iCNV.ref.ATAC[logi,])
  #add q7 averages to seurat object for each cell
  #average for observation cells
  q7.obs.av = colMeans(iCNV.obs.ATAC[logi,])
  q7.av = sapply(colnames(oligo.filt@assays$SCT),
                 function(y) {
                   if(y %in% names(c(q7.obs.av,q7.ref.av))){
                     c(q7.obs.av,q7.ref.av)[which(names(c(q7.obs.av,q7.ref.av))==y)]
                   }else{
                     NA
                   }
                 })
  names(q7.av) = colnames(oligo.filt@assays$SCT)
  oligo.filt$q7.av.ATAC = q7.av
  #find doublet score column to add to data frame for figure later
  md.doub = grepl("pANN_",colnames(oligo.filt@meta.data),fixed=T)
  if(grepl("oligo",samp)){
    p6 = FeaturePlot(oligo.filt, features = 'p1.av.ATAC', reduction = 'umap.rna', pt.size = 0.005) +
      scale_colour_gradientn(colours = pal, na.value = "grey90") +
      scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
    #boxplot showing 1p per cluster
    df = data.frame(cluster = oligo.filt@meta.data[,md.col], p1 = oligo.filt$p1.av.ATAC,
                    cell.type = oligo.filt@meta.data[,ct.col], doub.score = oligo.filt@meta.data[,md.doub])
    p7 = ggplot(df, aes(x=cluster, y=p1, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
      geom_hline(yintercept=(1-icnv.diff), color = "blue", linewidth=1.5)
  }else{
    p6 = FeaturePlot(oligo.filt, features = 'q7.av.ATAC', reduction = 'umap.rna', pt.size = 0.005) +
      scale_colour_gradientn(colours = pal, na.value = "grey90") +
      scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
    df = data.frame(cluster = oligo.filt@meta.data[,md.col], q7 = oligo.filt$q7.av.ATAC,
                    cell.type = oligo.filt@meta.data[,ct.col], doub.score = oligo.filt@meta.data[,md.doub])
    p7 = ggplot(df, aes(x=cluster, y=q7, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
      geom_hline(yintercept=(1+icnv.diff), color = "blue", linewidth=1.5) +
      coord_cartesian(ylim = c(0.95,1.13))
  }
}else{
  p6 = plot_spacer()
  p7 = plot_spacer()
}


###########################################
#### plot average of 19q genes in ATAC if oligo sample (since we see signal here but don't see in RNA)
###########################################

ref.pth = paste0(data.dir,"/analysis_",samp,"/infercnv_ATAC_",samp)

if(grepl("oligo",samp)){
  if(file.exists(paste0(ref.pth,"/infercnv.references.txt.gz"))){
    #input reference inferCNV
    iCNV.ref.ATAC = fread(paste0(ref.pth,"/infercnv.references.txt.gz"),data.table = F)
    rownames(iCNV.ref.ATAC) = iCNV.ref.ATAC[,1]
    iCNV.ref.ATAC = iCNV.ref.ATAC[,2:ncol(iCNV.ref.ATAC)]
    
    #input obs inferCNV
    iCNV.obs.ATAC = fread(paste0(ref.pth,"/infercnv.observations.txt.gz"),data.table = F)
    rownames(iCNV.obs.ATAC) = iCNV.obs.ATAC[,1]
    iCNV.obs.ATAC = iCNV.obs.ATAC[,2:ncol(iCNV.obs.ATAC)]
    
    #input gene order file
    gene.order.ATAC = fread(paste0(ref.pth,"/",samp,"_ATAC_cnv_ref.txt.gz"),data.table = F)
    
    q19.rows = which(pq.loc[,1]=="chr19" & grepl("^q",pq.loc[,4]))
    q19.st = pq.loc[q19.rows[1],2] #start of 7q
    q19.en = pq.loc[q19.rows[length(q19.rows)],3] #end of 7q
    #get all genes in this region on 7q
    q19.genes = gene.order.ATAC[gene.order.ATAC[,2]=="19" & gene.order.ATAC[,3]>q19.st & gene.order.ATAC[,4]<q19.en,1]
    #average inferCNV for 7q genes in reference cells
    logi = rownames(iCNV.ref.ATAC) %in% q19.genes
    #find number of 7q genes
    n7q = sum(logi)
    inds7q = which(logi)
    q19.ref.av = colMeans(iCNV.ref.ATAC[logi,])
    
    #add q19 averages to seurat object for each cell
    #average for observation cells
    q19.obs.av = colMeans(iCNV.obs.ATAC[logi,])
    q19.av = sapply(colnames(oligo.filt@assays$SCT),
                   function(y) {
                     if(y %in% names(c(q19.obs.av,q19.ref.av))){
                       c(q19.obs.av,q19.ref.av)[which(names(c(q19.obs.av,q19.ref.av))==y)]
                     }else{
                       NA
                     }
                   })
    names(q19.av) = colnames(oligo.filt@assays$SCT)
    oligo.filt$q19.av.ATAC = q19.av
    
    #find doublet score column to add to data frame for figure later
    md.doub = grepl("pANN_",colnames(oligo.filt@meta.data),fixed=T)
    p8 = FeaturePlot(oligo.filt, features = 'q19.av.ATAC', reduction = 'umap.rna', pt.size = 0.005) +
      scale_colour_gradientn(colours = pal, na.value = "grey90") +
      scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
    
    #boxplot showing 1p per cluster
    df = data.frame(cluster = oligo.filt@meta.data[,md.col], q19 = oligo.filt$q19.av.ATAC,
                    cell.type = oligo.filt@meta.data[,ct.col], doub.score = oligo.filt@meta.data[,md.doub])
    p9 = ggplot(df, aes(x=cluster, y=q19, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
      geom_hline(yintercept=(1-icnv.diff), color = "blue", linewidth=1.5)
  }else{
    p8 = plot_spacer()
    p9 = plot_spacer()
  }
}




####################################
####### SNV analysis from cb sniffer
####################################

#read in output file with cell barcodes and detected variants in RNA and ATAC if assay exists

if(file.exists(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                      "/", samp, "_RNA_counts_CB.tsv")) &&
   countLines(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                     "/", samp, "_RNA_counts_CB.tsv"))>1){
  cb.sn = read.table(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                            "/", samp, "_RNA_counts_CB.tsv"), sep = "\t", header = T,
                     colClasses = c("character",rep("numeric",2),rep("character",5),rep("numeric",3)))
}else{
  cb.sn = data.frame(chrm=NULL,start=NULL,end=NULL,ref=NULL,alt=NULL,
                     type=NULL,gene=NULL,barcode=NULL,ref_count=NULL,
                     alt_count=NULL,total_CB=NULL)
}
if(file.exists(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                      "/", samp, "_ATAC_counts_CB.tsv")) &&
   countLines(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                     "/", samp, "_ATAC_counts_CB.tsv"))>1){
  cb.sn = rbind(cb.sn, read.table(paste0(data.dir, "/analysis_", samp, "/cb_sniffer_", samp,
                                         "/", samp, "_ATAC_counts_CB.tsv"), sep = "\t", header = T,
                                  colClasses = c("character",rep("numeric",2),rep("character",5),rep("numeric",3))))
}

mut.uni = unique(cb.sn[,c(1,2,7)])
indMut = unlist(sapply(seq(nrow(mut.uni)), function(x) {
  rel.rows = cb.sn[cb.sn$chrm==mut.uni[x,1] & cb.sn$start==mut.uni[x,2],]
  if(sum(rel.rows$alt_count)>0){
    return(x)
  }else{
    return(NULL)
  }
}))

#create vector "IDH.mut" to add to seurat object with mutation info
#either "not_expressed" (no read covering mutation), "non_mutant" (read covering mutation, but no mutation detetced) or name of mutation
IDH.mut = rep("not_expressed", length(oligo.filt$orig.ident))
#if there are some mutations detected
if(length(indMut)>0){
  non_mut_cells = c()
  mut_cells = list()
  for(i in 1:length(indMut)){
    #find cells with mutant or non_mutant reads
    rel.rows = cb.sn[cb.sn$chrm==mut.uni[indMut[i],1] & cb.sn$start==mut.uni[indMut[i],2],]
    mut_cells[[i]] = rel.rows$barcode[rel.rows$alt_count>0]
    non_mut_cells = c(non_mut_cells, rel.rows$barcode[rel.rows$alt_count==0])
  }
  non_mut_cells = unique(non_mut_cells)
  mut_cells[[i+1]] = non_mut_cells
  #name lists after mutation
  names(mut_cells) = c(sapply(indMut, function(x) {paste(mut.uni[x,3], mut.uni[x,2], sep = "_")}), "non_mutant")
  
  #add mutation info to seurat object
  # add mutation info to oligo object to plot in umap
  IDH.mut = sapply(names(oligo.filt$orig.ident), function(x) {
    if(x %in% cb.sn$barcode){
      #if mutation is detected
      if(length(which(sapply(mut_cells, function(y) {x %in% y})))>0){
        #find which mutation
        return(names(mut_cells)[which(sapply(mut_cells, function(y) {x %in% y}))[1]])
      }else{
        return("non_mutant")
      }
    }else{
      return("not_expressed")
    }
  }, USE.NAMES = F)
  names(IDH.mut) = names(oligo.filt$orig.ident)
  oligo.filt$IDH.mut = IDH.mut
  
  mutants = unique(oligo.filt$IDH.mut)
  mutants = mutants[mutants!="not_expressed"]
  to_highlight = lapply(mutants, function(x) {
    names(oligo.filt$IDH.mut[which(oligo.filt$IDH.mut==x)])
  })
  names(to_highlight) = mutants
}else{
  IDH.mut[sapply(names(oligo.filt$orig.ident), function(x) {
    x %in% unique(cb.sn$barcode)
  })] = "non_mutant"
  names(IDH.mut) = names(oligo.filt$orig.ident)
  oligo.filt$IDH.mut = IDH.mut
  #colour all non_mutant
  non_mutant = unique(cb.sn$barcode)
  to_highlight = list(non_mutant = non_mutant)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(length(to_highlight))
sizs = rep(1.5,length(to_highlight))
sizs[names(to_highlight)=="non_mutant"] = 0.004

#count number of mutations detected in sample
mut.uni = unique(cb.sn[,c(1,2)])
nMut = sum(sapply(seq(nrow(mut.uni)), function(x) {
  rel.rows = cb.sn[cb.sn$chrm==mut.uni[x,1] & cb.sn$start==mut.uni[x,2],]
  if(sum(rel.rows$alt_count)>0){
    return(1)
  }else{
    return(0)
  }
}))
print(paste0(nMut," mutations detected"))
oligo.filt$IDH.mut = factor(oligo.filt$IDH.mut, ordered=T,
                            levels = c("not_expressed","non_mutant",
                                       unique(oligo.filt$IDH.mut[!(oligo.filt$IDH.mut %in% c("not_expressed","non_mutant"))])))

#umpa showing mutations detected
p10 = DimPlot(oligo.filt, group.by = 'IDH.mut', reduction = 'umap.rna', pt.size = 0.005,
             cells.highlight = to_highlight, cols.highlight = cols, raster = FALSE,
             order = c(c("not_expressed","non_mutant"),
                       unique(oligo.filt$IDH.mut[!(oligo.filt$IDH.mut %in% c("not_expressed","non_mutant"))]))) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")





################################
## Compute iCNV malignancy score
################################

# iCNV.score = sapply(names(oligo.filt$orig.ident), function(y) {
#   x = paste(y,oligo.filt$orig.ident[1],sep="_")
#   if(x %in% colnames(iCNV.ref)){
#     sum(abs(iCNV.ref[,which(colnames(iCNV.ref)==x)]-1))
#   }else if(x %in% colnames(iCNV.obs)){
#     sum(abs(iCNV.obs[,which(colnames(iCNV.obs)==x)]-1))
#   }else{
#     NA
#   }
# })
#in the end this simple (mean) count seems to work best
iCNV.score = sapply(names(oligo.filt$orig.ident), function(y) {
  x = paste(y,oligo.filt$orig.ident[1],sep="_")
  if(x %in% colnames(iCNV.ref)){
    sum(iCNV.ref[,which(colnames(iCNV.ref)==x)]<0.95 | iCNV.ref[,which(colnames(iCNV.ref)==x)]>1.05)/nrow(iCNV.ref)
  }else if(x %in% colnames(iCNV.obs)){
    sum(iCNV.obs[,which(colnames(iCNV.obs)==x)]<0.95 | iCNV.obs[,which(colnames(iCNV.obs)==x)]>1.05)/nrow(iCNV.obs)
  }else{
    NA
  }
})
df = cbind(df, iCNV.score)
names(df) = c(names(df)[1:(length(names(df))-1)],"iCNV.score")

#plot boxplot of CNV score by cluster
p11 = ggplot(df, aes(x=cluster, y=iCNV.score, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
  coord_cartesian(ylim = c(0,0.5))

#also doublet score
p12 = ggplot(df, aes(x=cluster, y=doub.score, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0))


#define malignant cells and plot
oligo.filt = define_malignant_cells(oligo.filt, cl.res, icnv.diff, sr.ref, pq.loc, iCNV.ref, iCNV.obs, data.dir)

#plot final umap showing which cells have been defined as malignant
p13 = DimPlot(oligo.filt, group.by = 'iCNV.malig.ident', reduction = 'umap.rna', pt.size = 0.005) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")


if(grepl("astro",samp)){
  print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p10 + p11 + p12 + p13 + plot_layout(ncol = 4))
  ggsave(paste0(data.dir, "/figures/umaps_malignant_summary_", samp, ".png"), 
         width=24, height=13)
}else{
  #oligos have extra figures of 19q in ATAC
  print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + plot_layout(ncol = 4))
  ggsave(paste0(data.dir, "/figures/umaps_malignant_summary_", samp, ".png"), 
         width=24, height=17)
}


#save object
saveRDS(oligo.filt, file = paste0(data.dir, "/analysis_", samp, "/oligo_filt_", samp, ".rds"))






