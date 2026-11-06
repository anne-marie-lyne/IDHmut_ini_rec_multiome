# code to cluster inferCNV ATAC output to identify CNA-based clones in each sample
#based on code from https://github.com/kokonech/mbOncoAberrations/blob/main/snATACseq/processAtacTumorSampleSV4.R

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp = as.character(args[2])


library(Seurat)
library(Signac)
library(infercnv, lib.loc="/mnt/beegfs/home/alyne/persistent/R_libs")
library(data.table)
library(GenomicRanges)
library(readxl)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#find subtype
cancer = unlist(strsplit(samp,"_"))[1]


#read in infercnv output
out.dir = paste0("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/infercnv_ATAC_reference_",samp)
# #read in seurat object for this sample
oligo = readRDS(paste0("/mnt/beegfs/home/alyne/stagein/data_tmp/numbat_inputs/numbat_",samp,"/oligo_filt_final_",samp,".rds"))
#read in gene order file from inferCNV analysis
gene.ord = fread(paste0("/mnt/beegfs/home/alyne/stagein/data_tmp/numbat_inputs/infercnv_atac_",cancer,
                        "/gene_order_",cancer,".tsv"), header=F, data.table = F, fill=TRUE)


obj <- readRDS(paste0(out.dir,"/run.final.infercnv_obj"))
obs = obj@expr.data
to.plot = as.matrix(t(obs)[unlist(obj@observation_grouped_cell_indices),])


#read in spreadsheet detailing CNA in each sample
CNVs <- as.data.frame(read_excel("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/summary_main_CNA_all_samples.xlsx", 
                                 col_types = "text", na = "NA"))


cell.states = c("astro-like","intermediate","neuro-like","OPC-like")
col.states = gg_color_hue(length(cell.states))
names(col.states) = cell.states


#if there are clonal CNV
if(!is.na(CNVs$CNA_use_clustering[which(CNVs$sample==samp)])){
  #find which chromosomes to use
  chr.to.plot = paste0("chr",unlist(strsplit(CNVs$CNA_use_clustering[which(CNVs$sample==samp)],";")))
  
  inds = match(colnames(to.plot), gene.ord$V1)
  col_split = gene.ord$V2[inds]
  
  #take only relevant chomrosomes
  to.plot.filt = to.plot[,col_split %in% chr.to.plot]
  col_split.filt = col_split[col_split %in% chr.to.plot]
  
  if(samp=="oligo_P0031_ini_S0060"){
    inds7 = which(col_split=="chr7")
    inds11 = which(col_split=="chr11")
    to.plot.filt = to.plot[,c(inds7[ceiling(length(inds7)/2):length(inds7)], inds11[ceiling(length(inds11)/2):length(inds11)])]
    col_split.filt = col_split[c(inds7[ceiling(length(inds7)/2):length(inds7)], inds11[ceiling(length(inds11)/2):length(inds11)])]
  }
  
  #make seurat object from icnv output
  assay.v5 <- CreateAssay5Object(data = t(to.plot.filt))
  # create a Seurat object based on this assay
  icnv <- CreateSeuratObject(assay.v5)
  # run standard Seurat pipeline to cluster inferCNV output
  all.genes <- rownames(icnv)
  icnv <- ScaleData(icnv, features = all.genes)
  icnv = RunPCA(icnv, features = all.genes)
  icnv <- FindNeighbors(icnv, dims = 1:10)
  icnv <- FindClusters(icnv, resolution = 0.1)
  if(samp %in% c("oligo_P0007_ini_S0015","astro_P0026_ini_S0050","oligo_P0031_ini_S0060")){
    icnv <- FindClusters(icnv, resolution = 0.2)
  }
  if(samp %in% c("oligo_P0005_rec1_S0012")){
    icnv <- FindClusters(icnv, resolution = 0.3)
  }
  icnv <- RunUMAP(icnv, dims = 1:10)
  ind = match(colnames(icnv),colnames(oligo))
  icnv$ICA_score_states = oligo$ICA_score_states[ind]
  #plot UMAPs coloured by clusters and lineage states
  DimPlot(icnv,reduction="umap",group.by = c("seurat_clusters","ICA_score_states")) + labs(x = "UMAP1",y = "UMAP2")
  ggsave(paste0(out.dir,"/infercnv_ATAC_ref_new_umap_",samp,".pdf"), width=10, height=4)
  
  
  #make figure showing proportions of cell states per clone
  df.cell.states = data.frame(icnv@meta.data[,c("seurat_clusters","ICA_score_states")]) %>%
    count(seurat_clusters, ICA_score_states)
  ggplot(df.cell.states, aes(x=seurat_clusters, y = n, fill = ICA_score_states)) + geom_bar(position="fill", stat="identity") +
    labs(fill='Cell state',x = "clone", y = "Proportion") + scale_fill_manual(values = col.states) +
    theme(text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0(out.dir,"/infercnv_ATAC_ref_new_clones_",samp,".pdf"), width=4, height=4)
  
  
  #make heatmap of clonal CNVs, split by clone
  df.row.anno = icnv@meta.data[,c("seurat_clusters","ICA_score_states")]
  df.row.anno = df.row.anno[order(df.row.anno$"seurat_clusters",df.row.anno$"ICA_score_states"),]
  
  inds = match(rownames(df.row.anno),rownames(to.plot.filt))
  to.plot.filt = to.plot.filt[inds,]
  
  row_split = df.row.anno$seurat_clusters
  
  col.cluster = rep("white", length(unique(df.row.anno$seurat_clusters)))
  names(col.cluster) = unique(df.row.anno$seurat_clusters)
  row_anno = rowAnnotation(df = df.row.anno, col = list(ICA_score_states = col.states, seurat_clusters = col.cluster),show_legend = c(F,T),
                           annotation_legend_param = list(ICA_score_states = list(title_gp = grid::gpar(fontsize = 18), labels_gp = grid::gpar(fontsize = 15)),
                                                          seurat_clusters = list(title_gp = grid::gpar(fontsize = 18), labels_gp = grid::gpar(fontsize = 15))))
  col_fun <- colorRamp2(c(0.92, 1, 1.08), c("blue", "white", "red"))
  ht <- Heatmap(
    to.plot.filt,
    col = col_fun,
    name = "iCNV",
    row_split = row_split,
    row_gap = unit(3, "mm"),
    cluster_rows = F,
    cluster_columns = F,
    border=T,
    column_split = factor(col_split.filt,levels = chr.to.plot),       # vector of cluster IDs
    column_title = chr.to.plot,
    column_title_gp = gpar(fontsize = 14),
    column_gap = unit(1, "mm"),
    left_annotation = row_anno,
    show_row_names = F,
    show_column_names = F,
    row_title = NULL,
    row_dend_width = unit(10, "cm"),  
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 18),labels_gp = grid::gpar(fontsize = 15))
  )
  pdf(paste0(out.dir,"/infercnv_ATAC_ref_new_by_clone_",samp,".pdf"), width = 10, height = 4)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}




