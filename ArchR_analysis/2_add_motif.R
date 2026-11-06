library(ArchR)
library(parallel)
library(data.table)
library(readxl)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(circlize)

#set subtype
cancer = "astro"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)

#ArchR parameters
set.seed(1)
addArchRLocking(locking = T)
addArchRThreads(threads = 6)
addArchRGenome("hg38")


#read in sample list
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]

oligo.ArchR = loadArchRProject(path = paste0(out.dir,"/project1"), force = T, showLogo = F)


######specific code to add annotations for jaspar2024
#already downloaded as no internet access from our cluster nodes
motif_mat = readRDS("/mnt/beegfs/home/alyne/persistent/motif_mat_JASPAR2024_human.rds")


#change names so that we have the TF name instead of the ID
adjust_labels <- function(vec) {
  dup_vals <- duplicated(vec) | duplicated(vec, fromLast = TRUE)
  vec[dup_vals] <- paste0(vec[dup_vals], "_", names(vec)[dup_vals])
  return(vec)
}
name_matrices <- sapply(motif_mat, function(pfm) TFBSTools::name(pfm))
deduplicated_name_matrices <- adjust_labels(name_matrices)
names(motif_mat) <- deduplicated_name_matrices[names(motif_mat)]
oligo.ArchR  <- addMotifAnnotations(oligo.ArchR, motifPWMs = motif_mat, annoName = "Motif", force = TRUE)

#save 
oligo.ArchR <- saveArchRProject(ArchRProj = oligo.ArchR, outputDirectory = paste0(out.dir,"/project1"), load = TRUE)


#add background peaks
oligo.ArchR <- addBgdPeaks(oligo.ArchR,force = TRUE)

#compute motif accessibilities using chromVAR
oligo.ArchR <- addDeviationsMatrix(
  ArchRProj = oligo.ArchR,
  peakAnnotation = "Motif",
  force = TRUE,
  threads = 6
)

#save
oligo.ArchR <- saveArchRProject(ArchRProj = oligo.ArchR, outputDirectory = paste0(out.dir,"/project1"), load = TRUE)




#get marker peaks for each malignant cell state
idxCells <- BiocGenerics::which(oligo.ArchR$cell_type %in% c("astro-like","OPC-like","intermediate","neuro-like"))
cellsSample <- oligo.ArchR$cellNames[idxCells]
oligo.ArchR.filt = oligo.ArchR[cellsSample, ]

peakAnno <- getPeakSet(oligo.ArchR.filt)

markerPeaks <- getMarkerFeatures(
  ArchRProj = oligo.ArchR.filt, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#filter for significant
heatmapPeaks <- markerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)


#heatmaps in Figure 4B
plotPDF(heatmapPeaks, name = paste0("Peak-Marker-Heatmap-",cancer), width = 8, height = 6,addDOC = FALSE)



