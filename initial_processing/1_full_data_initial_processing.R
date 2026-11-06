#script to run initial processing for single nucleus multiome data
#each section needs to be run independently per sample
#final output is filtered Seurat object per sample


rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
pbs.arrayid = as.numeric(args[3])


setwd(data.dir)

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(DoubletFinder)
library(RColorBrewer)
library(readxl)
library(future)
library(destiny)
options(future.globals.maxSize = 100 * 1024 ^ 3)


#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric"), 
                     stringsAsFactors = F, fill = T)
samples = sample.info[,1]
samples.newnames = sample.info[,2]
samples.new.sh = sapply(strsplit(samples.newnames,"_"),
                        function(x) {paste(x[2:4],collapse = "_")})
sample.date = sample.info[,3]


##########################
#### load in raw data and make Seurat objects
##########################

library(hdf5r)

#the 10x hdf5 file contains both data types.
inputdata.10x <- Read10X_h5(paste0("analysis_",samples.newnames[pbs.arrayid],"/CR_output/",
                                   "filtered_feature_bc_matrix.h5"))

# extract RNA and ATAC data if both present, otherwise just RNA
if(length(inputdata.10x)==2){
  rna_counts <- inputdata.10x$`Gene Expression`
}else{
  rna_counts <- inputdata.10x
}

# Create Seurat object
oligo <- CreateSeuratObject(counts = rna_counts, min.cells = 3)
oligo[["percent.mt"]] <- PercentageFeatureSet(oligo, pattern = "^MT-")
oligo =  PercentageFeatureSet(oligo, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                              col.name = "percent.ribo")

if(length(inputdata.10x)==2){
  atac_counts <- inputdata.10x$Peaks
  # Now add in the ATAC-seq data
  # only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

  frag.file <- paste0("analysis_",samples.newnames[pbs.arrayid],"/CR_output/","atac_fragments.tsv.gz")

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    # genome = 'hg38',
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
  )
  oligo[["ATAC"]] <- chrom_assay

  DefaultAssay(oligo) <- "ATAC"
  oligo <- NucleosomeSignal(oligo)
  oligo <- TSSEnrichment(oligo)
}else{
  print("Only GEX")
}

#add sample name as orig.ident
oligo$orig.ident = rep(samples.newnames[pbs.arrayid], length(oligo@meta.data$orig.ident))
oligo@active.ident = as.factor(oligo$orig.ident)

print(paste0("sample ",samples.newnames[pbs.arrayid]," has ",length(oligo@meta.data$orig.ident)," cells"))

#save before any filtering
saveRDS(oligo, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                      "/oligo_seurat_", samples.newnames[pbs.arrayid], ".rds"))



#######################################
## Now separately run code to make pre-filtering QC plots (make_figures_QC.R)
#######################################



##############################
## Do QC filtering and doublet detection
##############################


oligo = readRDS(file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                      "/oligo_seurat_", samples.newnames[pbs.arrayid], ".rds"))

oligo[["RNA"]] <- as(object = oligo[["RNA"]], Class = "Assay")
N_prefilt = ncol(oligo)

#get 10x multiplet rate
tab.multiplet = as.data.frame(read_excel(paste0(data.dir,"/multiplet_rate_table.xlsx"),
                                         col_types = "numeric"))
print(head(tab.multiplet,11))
doub.rate.10x = tab.multiplet[which.min(abs(tab.multiplet[,3]-N_prefilt)),1]
print(doub.rate.10x)

#check if there are already any metadata columns from doublet detection and remove
ind1 = which(grepl("DF.classifications_",colnames(oligo@meta.data),fixed=T))
if(length(ind1)>0){
  oligo@meta.data = oligo@meta.data[,-ind1]
}
ind2 = which(grepl("pANN_",colnames(oligo@meta.data),fixed=T))
if(length(ind2)>0){
  oligo@meta.data = oligo@meta.data[,-ind2]
}

#first filter on MT%
mt.thresh = 5
oligo <- subset(oligo, subset = percent.mt < mt.thresh)

#also filters on other quality metrics depending on which assays present
min.count.ATAC = 1000
min.count.RNA = 1000
max.count.RNA = 100000
min.feat.RNA = 500

if("nCount_ATAC" %in% colnames(oligo@meta.data)){
  oligo$log10_nCount_ATAC = log10(oligo$nCount_ATAC+1)
  oligo <- subset(oligo, subset = nCount_RNA > min.count.RNA &
                    nCount_RNA < max.count.RNA &
                    nFeature_RNA > min.feat.RNA &
                    TSS.enrichment > 2 &
                    nCount_ATAC > min.count.ATAC)
}else{
  oligo <- subset(oligo, subset = nCount_RNA > min.count.RNA &
                    nCount_RNA < max.count.RNA &
                    nFeature_RNA > min.feat.RNA)
}


# run initial SCT processing, don't regress anything in order to match doublet finder
oligo <- SCTransform(oligo, verbose = FALSE)
oligo <- RunPCA(oligo)
oligo <- FindNeighbors(oligo, dims = 1:30, verbose = FALSE)
oligo <- FindClusters(oligo, verbose = FALSE, graph.name = "SCT_snn")

#count number of cells before doublet filtering and compute expected number of doublets
N_postfilt = ncol(oligo)
print(N_postfilt)
Nd.10x = ceiling(N_postfilt*doub.rate.10x/100)
print(Nd.10x)


#DoubletFinder
sweep.res.list_oligo <- paramSweep(oligo, PCs = 1:30, sct = TRUE)
sweep.stats_oligo <- summarizeSweep(sweep.res.list_oligo, GT = FALSE)
bcmvn_oligo <- find.pK(sweep.stats_oligo)

#find pK for this dataset
pK.opt = bcmvn_oligo$pK[which.max(bcmvn_oligo$BCmetric)]
pK.opt = as.numeric(levels(pK.opt))[pK.opt]


## Run DoubletFinder
oligo <- doubletFinder(oligo, PCs = 1:30, pN = 0.25, pK = pK.opt,
                          nExp = Nd.10x, reuse.pANN = FALSE, sct = TRUE)


#run SCT processing regressing %MT and Ncount
oligo <- SCTransform(oligo, verbose = FALSE, vars.to.regress = c('percent.mt','nCount_RNA'))
oligo <- RunPCA(oligo)
oligo <- FindNeighbors(oligo, dims = 1:30, verbose = FALSE)
oligo <- FindClusters(oligo, verbose = FALSE, graph.name = "SCT_snn")
oligo <- RunUMAP(oligo, dims = 1:30, reduction.key = "UMAP_", reduction.name = "umap.rna")


#plot umaps of doublet finder scores to see if there are any remaining
#clusters which are likely to be doublets
ind1 = grepl("DF.classifications_",colnames(oligo@meta.data),fixed=T)
ind2 = grepl("pANN_",colnames(oligo@meta.data),fixed=T)


#get threshold using 10x multiplet rate
thresh.10x = sort(oligo@meta.data[,ind2],decreasing = T)[Nd.10x]


df = data.frame(Cluster = oligo$SCT_snn_res.0.8, Doublet.score = oligo@meta.data[,ind2],
                Class = oligo@meta.data[,ind1], N.feature = oligo$nFeature_RNA)
df = rbind(df, data.frame(Cluster = "All", Doublet.score = oligo@meta.data[,ind2],
                          Class = oligo@meta.data[,ind1], N.feature = oligo$nFeature_RNA))

p1 = DimPlot(oligo, group.by = "SCT_snn_res.0.8", reduction = 'umap.rna') +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
p2 = FeaturePlot(oligo, features = colnames(oligo@meta.data)[ind2], reduction = 'umap.rna',
                 pt.size = 0.005) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
p3 = DimPlot(oligo, group.by = colnames(oligo@meta.data)[ind1], reduction = 'umap.rna',
                 pt.size = 0.005) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
p4 = ggplot(df, aes(x=Cluster, y=Doublet.score)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept=thresh.10x, color = "green", linewidth=1.5)
p5 <- ggplot(df, aes(x=Doublet.score)) + geom_density() +
  geom_vline(xintercept=thresh.10x, color = "green", linewidth=1.5) +
  ggtitle(paste0(N_postfilt," cells, 10x doublet rate ",doub.rate.10x))
p6 = ggplot(df, aes(x=N.feature, y=Doublet.score)) + geom_point()
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3))
ggsave(paste0("figures/umaps_doublet_prefilt_", samples.newnames[pbs.arrayid], ".pdf"),
       width = 15, height = 9)


#find which cells to keep from doubletfinder classification
cells.to.keep = names(oligo$orig.ident)[which(oligo@meta.data[,ind1]=="Singlet")]

#filter doublets
oligo = subset(oligo, cells=cells.to.keep)

# re-run initial processing, regressing
oligo <- SCTransform(oligo, verbose = FALSE, vars.to.regress = c('percent.mt', 'nCount_RNA'))
oligo <- RunPCA(oligo)
oligo <- FindNeighbors(oligo, dims = 1:30, verbose = FALSE)
oligo <- FindClusters(oligo, verbose = FALSE, graph.name = "SCT_snn")
oligo <- FindClusters(oligo, verbose = FALSE, graph.name = "SCT_snn", resolution = 1.2)
oligo <- FindClusters(oligo, verbose = FALSE, graph.name = "SCT_snn", resolution = 0.3)
oligo <- RunUMAP(oligo, dims = 1:30, reduction.key = "UMAP_", reduction.name = "umap.rna")

#save output
saveRDS(oligo, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                             "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))


#plot umaps again of doublet finder scores to see if there are any remaining
#clusters which are likely to be doublets
ind2 = grepl("pANN_",colnames(oligo@meta.data),fixed=T)
p1 = FeaturePlot(oligo, features = colnames(oligo@meta.data)[ind2], reduction = 'umap.rna', pt.size = 0.005) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
p2 = DimPlot(oligo, group.by = "SCT_snn_res.1.2", reduction = 'umap.rna') +
  scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
print(p1 | p2)
ggsave(paste0("figures/umaps_doublet_postfilt_", samples.newnames[pbs.arrayid], ".pdf"),
       width = 12, height = 5)





################################
## dim reduc, clustering etc on ATAC assay using signac and get joint projection 
##(don't really use in the final paper...)
################################

oligo = readRDS(file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                              "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))

if("ATAC" %in% names(oligo@assays)){

  DefaultAssay(oligo) <- "ATAC"
  oligo <- RunTFIDF(oligo)
  oligo <- FindTopFeatures(oligo, min.cutoff = 'q0')
  oligo <- RunSVD(oligo)
  oligo <- FindNeighbors(object = oligo, reduction = 'lsi', dims = 2:50)
  oligo <- FindClusters(object = oligo, verbose = FALSE, algorithm = 3, graph.name = "ATAC_snn")
  oligo <- RunUMAP(oligo, reduction = 'lsi', dims = 2:50,
                   reduction.name = "umap.atac", reduction.key = "atacUMAP_")

  #compute WNN graph representing weighted combination of RNA and ATAC
  oligo <- FindMultiModalNeighbors(oligo, reduction.list = list("pca", "lsi"),
                                   dims.list = list(1:50, 2:50))
  oligo <- FindClusters(oligo, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  oligo <- RunUMAP(oligo, nn.name = "weighted.nn", reduction.name = "wnn.umap",reduction.key = "wnnUMAP_")
}


if("ATAC" %in% names(oligo@assays)){
  #umaps of rna atac and combined
  p1 = DimPlot(oligo, group.by = "SCT_snn_res.0.8", reduction = 'umap.rna') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p2 = DimPlot(oligo, group.by = "ATAC_snn_res.0.8", reduction = 'umap.atac') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p3 = DimPlot(oligo, group.by = "wsnn_res.0.8", reduction = 'wnn.umap') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  print(p1 | p2 | p3)
  ggsave(paste0("figures/umaps_", samples.newnames[pbs.arrayid], ".pdf"),
         width = 12, height = 4)

  p1 = DimPlot(oligo, group.by = "SCT_snn_res.0.8", reduction = 'umap.rna') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p2 = DimPlot(oligo, group.by = "ATAC_snn_res.0.8", reduction = 'umap.rna') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p3 = DimPlot(oligo, group.by = "SCT_snn_res.0.8", reduction = 'umap.atac') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  p4 = DimPlot(oligo, group.by = "ATAC_snn_res.0.8", reduction = 'umap.atac') +
    scale_x_continuous(name = "UMAP 1") + scale_y_continuous(name = "UMAP 2")
  print((p1 | p2) / (p3 | p4))
  ggsave(paste0("figures/umaps_RNA_ATAC_", samples.newnames[pbs.arrayid], ".pdf"), width = 9, height = 8)
}

#make list of cells which passed filters
fwrite(list(names(oligo$orig.ident)),
            file = paste0(data.dir,"/analysis_",samples.newnames[pbs.arrayid],
                          "/cells_pass_filt_both_",samples.newnames[pbs.arrayid],".txt.gz"),
            row.names = F, col.names = F, quote = F)

#save output
saveRDS(oligo, file = paste0(data.dir, "/analysis_", samples.newnames[pbs.arrayid],
                             "/oligo_filt_", samples.newnames[pbs.arrayid], ".rds"))







