library(ArchR)
library(parallel)
library(data.table)
library(readxl)

#set cancer subtype
cancer = "astro"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)

set.seed(1)

#set ArchR parameters
addArchRLocking(locking = T)
addArchRThreads(threads = 6)
addArchRGenome("hg38")


#read in sample list
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]


#read in cells for each sample (which passed scenicplus QC)
inputFiles = paste0("stagein/scenicplus_",cancer,"/analysis_",samps,"/CR_output/atac_fragments.tsv.gz")


#read in rss file to get cell names which have already passed filtering by scenic plus
cells.sp = fread(paste0("stageout/scenicplus_op_newdb_filt_",cancer,"/figures_scenicplus_newdb_",
                        cancer,"/direct_gene_based_AUC_",cancer,".tsv.gz"), data.table=F, header=T)[,1]
cells.list = lapply(samps, function(x) {
  cells.rel = cells.sp[grepl(x, cells.sp)]
  sapply(cells.rel, function(y) {paste(unlist(strsplit(y,"-"))[1:2], collapse="-")})
})
names(cells.list) = samps



#create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = samps,
  outputNames = paste0(out.dir,"/",samps),
  minTSS = 0, #set these to zero as filtering has already been done
  minFrags = 0, 
  addTileMat = T,
  addGeneScoreMat = T,
  validBarcodes = cells.list,
  force = F
)


Nsamps.use = length(samps)
ArrowFiles = paste0(out.dir,"/",samps,".arrow")


#create archr project
oligo.ArchR <- ArchRProject(ArrowFiles = ArrowFiles, 
  outputDirectory = paste0(out.dir,"/project1"), copyArrows = F)


#compute LSI dimensino reduction
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

#clustering
oligo.ArchR <- addClusters(
  input = oligo.ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2, force = T
)

#UMAP
oligo.ArchR <- addUMAP(
  ArchRProj = oligo.ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


#add peaks from pycistopic so have one consistent set for all 
peaks <- import(paste0("stageout/pycisTopic_output_",cancer,"/pycistopic_peaks_",cancer,".bed"))
oligo.ArchR <- addPeakSet(ArchRProj = oligo.ArchR, peakSet = peaks, force = T)
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


#read in cell annotations
cell.annot = fread(paste0("stagein/scenicplus_",cancer,"/cell_annotations.tsv.gz"),header = T, data.table = F)
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
oligo.ArchR$treatment = sapply(oligo.ArchR$Sample, function(x) {
  TMZ = samp.md$TMZ[which(samp.md$Sample==x)]
  PCV = samp.md$PCV[which(samp.md$Sample==x)]
  if(TMZ==0 && PCV==0){
    return("none")
  }else{
    if(TMZ==1){
      return("TMZ")
    }else{
      return("PCV")
    }
  }
})


#save
oligo.ArchR <- saveArchRProject(ArchRProj = oligo.ArchR, outputDirectory = paste0(out.dir,"/project1"), load = T)
  

# oligo.ArchR <- loadArchRProject(path = paste0(out.dir,"/project1"), force = T, showLogo = TRUE)

#integrate using harmony
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

#ggplot colour scale function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#make colour palette matching RNA figures
cell.types = c("Astrocyte","Endothelial","GABAergic","Glutamatergic","Lymphoid","Malignant","Myeloid", "Oligodendrocyte","OPC")
col.pal = gg_color_hue(length(cell.types))
names(col.pal) = cell.types

#add metadata with correct cell type names
oligo.ArchR$cell_type_match = oligo.ArchR$cell_type
oligo.ArchR$cell_type_match[which(oligo.ArchR$cell_type_match=="Astro")] = "Astrocyte"
oligo.ArchR$cell_type_match[which(oligo.ArchR$cell_type_match=="Oligo")] = "Oligodendrocyte"
oligo.ArchR$cell_type_match[which(oligo.ArchR$cell_type_match=="Macrophage")] = "Myeloid"
oligo.ArchR$cell_type_match[which(oligo.ArchR$cell_type_match %in% c("astro-like","OPC-like","intermediate","neuro-like"))] = "Malignant"


#plot UMAPs in Figure supp figure 4A
p1 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony",labelMeans=F)
p2 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP",labelMeans=F)
p3 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "cell_type_match", embedding = "UMAPHarmony",labelMeans=F, pal = col.pal)
p4 = plotEmbedding(ArchRProj = oligo.ArchR, colorBy = "cellColData", name = "cell_type_match", embedding = "UMAP",labelMeans=F, pal = col.pal)
plotPDF(p1,p2,p3,p4, name = paste0("Plot-UMAP2Harmony-Sample-Clusters_",cancer,".pdf"), 
        ArchRProj = oligo.ArchR, addDOC = FALSE, width = 10, height = 10)

