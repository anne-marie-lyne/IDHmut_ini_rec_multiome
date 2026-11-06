library(ArchR)
library(parallel)
library(data.table)
library(readxl)

cancer = "astro"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)

set.seed(1)

addArchRLocking(locking = T)

addArchRThreads(threads = 6)

addArchRGenome("hg38")

#read in sample list
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]

#read in cells for each sample (which passed scenicplus QC)
inputFiles = paste0("stagein/scenicplus_",cancer,"/analysis_",samps,"/CR_output/atac_fragments.tsv.gz")

#read in rss file to get cell names which have already passed filtering
cells.sp = fread(paste0("stageout/scenicplus_op_newdb_filt_",cancer,"/figures_scenicplus_newdb_",
             cancer,"/direct_gene_based_AUC_",cancer,".tsv.gz"), data.table=F, header=T)[,1]
cells.list = lapply(samps, function(x) {
  cells.rel = cells.sp[grepl(x, cells.sp)]
  sapply(cells.rel, function(y) {paste(unlist(strsplit(y,"-"))[1:2], collapse="-")})
})
names(cells.list) = samps


#read in cell annotations
cell.annot = fread(paste0("stagein/scenicplus_",cancer,"/cell_annotations.tsv.gz"),header = T, data.table = F)

#add column to cell.annot so cell names match
cell.annot$names_scenic = paste(cell.annot$V1, cell.annot$sample, sep="___")


#get barcodes of malignant cells
malig.types = c("astro-like","intermediate","OPC-like","neuro-like")
cells.list.malig = lapply(cells.list, function(x) {
  scenic.names = names(x)
  samp = unlist(strsplit(scenic.names[1],"___"))[2]
  #cells from this sample which are malignant
  malig.cells = cell.annot$names_scenic[(cell.annot$cell_type %in% malig.types) & cell.annot$sample == samp]
  x[scenic.names %in% malig.cells]
})


#create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = samps,
  outputNames = paste0(out.dir,"/",samps,"_malig"),
  minTSS = 0, #set these to zero as filtering has already been done
  minFrags = 0,
  addTileMat = T,
  addGeneScoreMat = T,
  validBarcodes = cells.list.malig,
  force = T
)

Nsamps.use = length(samps)
ArrowFiles = paste0(out.dir,"/",samps,"_malig.arrow")

#create archr project
oligo.ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0(out.dir,"/project3"),
  copyArrows = F
)

oligo.ArchR <- addIterativeLSI(
  ArchRProj = oligo.ArchR,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

oligo.ArchR <- addClusters(
  input = oligo.ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2, force = T
)

oligo.ArchR <- addUMAP(
  ArchRProj = oligo.ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

#add peaks from pycistopic
peaks <- import(paste0("stageout/pycisTopic_output_",cancer,"/pycistopic_peaks_",cancer,".bed"))
oligo.ArchR <- addPeakSet(
  ArchRProj = oligo.ArchR,
  peakSet = peaks, force = T
)
oligo.ArchR <- addPeakMatrix(oligo.ArchR)


#add meta data
oligo.ArchR$Sample.number = sapply(oligo.ArchR$Sample, function(x) {
  unlist(strsplit(x,"_"))[4]
})
oligo.ArchR$Patient.number = sapply(oligo.ArchR$Sample, function(x) {
  unlist(strsplit(x,"_"))[2]
})
oligo.ArchR$ini_rec = sapply(oligo.ArchR$Sample, function(x) {
  unlist(strsplit(x,"_"))[3]
})


#swap around cell names to match ArchR
cell.names.ArchR = sapply(cell.annot$V1, function(x) {
  samp = unlist(strsplit(x, "-"))[3]
  cell.name = paste(unlist(strsplit(x, "-"))[1:2],collapse="-")
  paste0(samp,"#",cell.name)
})
oligo.ArchR$cell_type = sapply(oligo.ArchR$cellNames, function(x) {cell.annot$cell_type[which(cell.names.ArchR==x)]})


#read in grade and treatment info
samp.md = as.data.frame(read_excel("persistent/sample_summary_astro_oligo_anon.xlsx"))

oligo.ArchR$grade = sapply(oligo.ArchR$Sample, function(x) {
  samp.md$Grade[which(samp.md$Sample==x)]
})
oligo.ArchR$grade_comb = sapply(oligo.ArchR$grade, function(x) {
  if(x==4){
    return("4")
  }else{
    return("2_3")
  }
})
if(cancer=="astro"){
  oligo.ArchR$treatment = sapply(oligo.ArchR$Sample, function(x) {
    samp.md$TMZ[which(samp.md$Sample==x)]
  })
}else{
  oligo.ArchR$treatment = sapply(oligo.ArchR$Sample, function(x) {
    samp.md$PCV[which(samp.md$Sample==x)]
  })
}


oligo.ArchR <- saveArchRProject(ArchRProj = oligo.ArchR, outputDirectory = paste0(out.dir,"/project3"), load = T)
# oligo.ArchR <- loadArchRProject(path = paste0(out.dir,"/project3"), force = T, showLogo = TRUE)


#integrate with harmony
oligo.ArchR <- addHarmony(
  ArchRProj = oligo.ArchR,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

oligo.ArchR <- addUMAP(
  ArchRProj = oligo.ArchR, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


#ggplot2 palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#ake colour palette matching RNA figures
cell.types = c("astro-like","intermediate","neuro-like","OPC-like")
col.pal = gg_color_hue(length(cell.types))
names(col.pal) = cell.types


##UMAPs in Figure 4A
p1 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony",labelMeans=F)
p2 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP",labelMeans=F)
p3 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony",labelMeans=F, pal = col.pal)
p4 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP",labelMeans=F, pal = col.pal)
plotPDF(p1,p2,p3,p4, name = paste0("Plot-UMAP2Harmony-Sample-Clusters_",cancer,"_malignant.pdf"), 
        ArchRProj = oligo.ArchR, addDOC = FALSE, width = 10, height = 10)



