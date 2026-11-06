#make violin plots pre-filtering
#subsample to plot 1/4 of points

data.dir = "/data/users/alyne/gliomult/svn/analyse/results"
samp.file = "/data/users/alyne/scrnaseq/oligo_10x_multiome/oligo_astro_samples_all.txt"

setwd(data.dir)

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(RColorBrewer)

batch = 3 # divide samples into three batches so can be seen in figures
fig = 0 # set fig=0 for prefiltering, fig = 1 for postfiltering

#define which seurat objects and figure names to use
if(fig==0){
  seu.obj = "seurat"
  fig.str = "prefilter"
}else{
  seu.obj = "filt"
  fig.str = "postfilter"
}

#define which samples to include in each batch
st = c(1,23,41)
en = c(22,40,59)

#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric"), 
                         stringsAsFactors = F, fill = T)
samples.newnames = sample.info[,2]
samples.newnames = sort(samples.newnames)


for(batch in 1:3){
  print(paste0("batch ",batch))
  rel.samps = samples.newnames[st[batch]:en[batch]]
  oligo.list = list()
  for(i in 1:length(rel.samps)){
    print(i)
    oligo.list[[i]] = readRDS(paste0("analysis_", rel.samps[i],"/oligo_",seu.obj,"_", rel.samps[i], ".rds"))
  }
  df.list = lapply(oligo.list, function(x) {
    if("nCount_ATAC" %in% colnames(x@meta.data)){
      nCount_ATAC = x$nCount_ATAC
    }else{
      nCount_ATAC = NA
    }
    if("TSS.enrichment" %in% colnames(x@meta.data)){
      TSS.enrichment = x$TSS.enrichment
    }else{
      TSS.enrichment = NA
    }
    data.frame(samp = x$orig.ident,
               patient = unlist(strsplit(x$orig.ident[1],"_"))[2],
               type = unlist(strsplit(x$orig.ident[1],"_"))[3],
               nCount_ATAC = nCount_ATAC+1,
               TSS.enrichment = TSS.enrichment,
               nFeature_RNA = x$nFeature_RNA,
               nCount_RNA = x$nCount_RNA+1,
               percent.mt = x$percent.mt)
  })
  df1 = do.call("rbind", df.list)
  df2 = df1[sample(nrow(df1),floor(nrow(df1)/2),replace=F),]
  plot.list = list()
  j=1
  for(i in 4:ncol(df1)){
    p1 = ggplot(df2, aes(x=samp, y=.data[[colnames(df1)[i]]], fill=samp)) + geom_violin(trim=F) + 
      geom_point(shape='.', position = position_jitter(seed = 1, width = 0.2)) +
      theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = rel(3))) +
      theme(axis.text.y = element_text(size = rel(3)), axis.title.y = element_text(size = rel(3)))
    if(i==4 || i==7){
      p1 = p1 + scale_y_log10()
    }
    if(i==5 || i==8){
      p1 = p1 + coord_cartesian(ylim = c(0, 40))
    }
    plot.list[[j]] = p1
    j = j + 1
  }
  p1 = wrap_plots(plot.list, ncol=1)
  print(p1)
  ggsave(paste0("figures/QC_plot_",fig.str,"_",batch,".png"), width = 20, height = 50, limitsize = FALSE)
  
  #########################################################
  #if pre-filtering, replot figures after just MT filtering
  #########################################################
  for(mt.thresh in c(3,5)){
    if(fig==0){
      oligo.filt = lapply(oligo.list, function(x) {
        subset(x, subset = percent.mt < mt.thresh)
      })
      df.list = lapply(oligo.filt, function(x) {
        if("nCount_ATAC" %in% colnames(x@meta.data)){
          nCount_ATAC = x$nCount_ATAC
        }else{
          nCount_ATAC = NA
        }
        if("TSS.enrichment" %in% colnames(x@meta.data)){
          TSS.enrichment = x$TSS.enrichment
        }else{
          TSS.enrichment = NA
        }
        data.frame(samp = x$orig.ident,
                   patient = unlist(strsplit(x$orig.ident[1],"_"))[2],
                   type = unlist(strsplit(x$orig.ident[1],"_"))[3],
                   nCount_ATAC = nCount_ATAC+1,
                   TSS.enrichment = TSS.enrichment,
                   nFeature_RNA = x$nFeature_RNA,
                   nCount_RNA = x$nCount_RNA+1,
                   percent.mt = x$percent.mt)
      })
      df1 = do.call("rbind", df.list)
      df2 = df1[sample(nrow(df1),floor(nrow(df1)/2),replace=F),]
      plot.list = list()
      j = 1
      for(i in 4:ncol(df1)){
        p1 = ggplot(df2, aes(x=samp, y=.data[[colnames(df1)[i]]], fill=samp)) + geom_violin(trim=F) + 
          geom_point(shape='.', position = position_jitter(seed = 1, width = 0.2)) +
          theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_text(size = rel(3))) +
          theme(axis.text.y = element_text(size = rel(3)), axis.title.y = element_text(size = rel(3)))
        if(i==4){
          p1 = p1 + geom_hline(yintercept=500, linewidth=1, color="green")
        }
        if(i==5){
          p1 = p1 + geom_hline(yintercept=2, linewidth=1, color="green")
        }
        if(i==6){
          p1 = p1 + geom_hline(yintercept=500, linewidth=1, color="green")
        }
        if(i==4 || i==7){
          p1 = p1 + scale_y_log10()
        }
        if(i==5){
          p1 = p1 + coord_cartesian(ylim = c(0, 30))
        }
        if(i==8){
          p1 = p1 + coord_cartesian(ylim = c(0, 5))
        }
        plot.list[[j]] = p1
        j = j + 1
      }
      p1 = wrap_plots(plot.list, ncol=1)
      print(p1)
      ggsave(paste0("figures/QC_plot_MTfilt",mt.thresh,"_",batch,".png"), 
             width = 20, height = 50, limitsize = FALSE)
    }
  }
}






