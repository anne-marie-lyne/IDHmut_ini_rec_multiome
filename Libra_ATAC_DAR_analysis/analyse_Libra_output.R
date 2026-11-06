#script to analyse LIBRA output
#plot heatmaps of DARs and enriched motifs for figure 5
rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
folder = as.character(args[1])
str.ana = as.character(args[2])
cancer = as.character(args[3])
analysis_name = as.character(args[4])
archr.proj = as.character(args[5])
grp.forward = as.character(args[6])
grp.backward = as.character(args[7])


########################
####plot heatmap of DARs in Figure 5
########################


#choose one of the following and comment out the others for the comparisons presented in Fig 5

# compare grade 4 IDH-A TMZ-treated to untreated
folder = "diff_output_merged"
str.ana = "TMZ_gr4"
samp.md.col = "TMZ_gr4"
cancer = "astro"
analysis_name = "astrocytoma gr4 - treatment"
archr.proj = "project2"
grp.forward = "TMZ4"
grp.backward = "noTMZ4"

# compare grade 4 untreated IDH-A to untreated grades 2/3
folder = "diff_output_merged"
str.ana = "gr_noTr"
samp.md.col = "gr_noTr"
cancer = "astro"
analysis_name = "astrocytoma untreated - grade"
archr.proj = "project2"
grp.forward = "4"
grp.backward = "2_3"

# compare grade 3 IDH-O PCV-treated to untreated
folder = "diff_output_merged"
str.ana = "PCV_gr3"
samp.md.col = "PCV_gr3"
cancer = "oligo"
analysis_name = "oligodendroglioma gr3 - PCV"
archr.proj = "project2"
grp.forward = "PCV3"
grp.backward = "noPCV3"

if(folder == "diff_output_merged"){
  str.mer = "_merged"
}else{
  str.mer = ""
}

library(ArchR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)
library(readxl)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(pvclust)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(scico)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols.TF = gg_color_hue(2)

OUT_LIBRA = paste0("/mnt/beegfs/home/alyne/stageout/ArchR_output_",cancer,"/",folder)
setwd(OUT_LIBRA)
dir.create("figures_Libra")
CODE_PATH = "/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas"
# import results 
results_import <- readRDS(file.path(OUT_LIBRA, paste0(archr.proj,"_libra_diff_",str.ana,"_",
                                                      grp.backward,"_",grp.forward,"_No_control_DESeq2_Wald.rds")))

#read in grade and treatment info
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))

#functions scripts
source_all_r_files <- function(path_list){
  origin_wd = getwd()
  lapply(
    path_list, function(x) {
      setwd(x)
      file.sources = list.files(pattern =  "*.r$")
      sapply(file.sources,source,.GlobalEnv)
    }
  )
  setwd(origin_wd)
}
source_all_r_files(c(file.path(CODE_PATH, "Miscellaneous"), file.path(CODE_PATH, "Differential_analysis")))



# convert results to input for plotting function 
se_list_astro_both <- lapply(colnames(results_import$Peaks), function(col) {
  list(Peaks = results_import$Peaks[,col], MotifsUp = results_import$MotifsUp[,col], MotifsDown = results_import$MotifsDown[,col], 
       Genes = NULL)#results_genes_astro$results[results_genes_astro$results$cell_type == col,])
})
names(se_list_astro_both) <- colnames(results_import$Peaks)

motifs_up_astro_df  <- lapply(names(se_list_astro_both), function(name) getDFfromSESample(se_list_astro_both[[name]]$MotifsUp))
names(motifs_up_astro_df) <- names(se_list_astro_both)
motifs_down_astro_df  <- lapply(names(se_list_astro_both), function(name) getDFfromSESample(se_list_astro_both[[name]]$MotifsDown))
names(motifs_down_astro_df) <- names(se_list_astro_both)


#dendrograms, heatmaps and PCA
full_mat_scaled_list <- lapply(results_import$matrices, function(i) {
  if(!is.null(i)){
    a <- scale(t(i))
    print("done")
    return(a)
  }else{NULL}
})
names(full_mat_scaled_list) <- names(results_import$matrices)

cols.bin = c("1" = "grey30", "0" = "gray90")
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
# cols.ini_rec = c("ini" = "violet", "rec1" = "purple")
ann_colors = list(grade = cols.grade, ini_rec = cols.ini_rec, TMZ = cols.bin, PCV = cols.bin, RT = cols.bin)


cell.type = "malignant"

# Extract significant genes
sig_variables <- results_import$Libra %>% filter(cell_type == cell.type, p_val_adj < 0.05, abs(avg_logFC) > 0.5) %>% pull(gene)

if(length(sig_variables)>1){
  # Subset matrix for significant genes
  mat <- t(results_import$matrices[[cell.type]][sig_variables, ])
  log_mat <- log1p(mat)
  
  # Perform pvclust bootstrapped clustering on samples (columns)
  # Transpose is needed: rows = samples, columns = genes
  pv_result <- pvclust(t(log_mat), method.hclust = "ward.D2", method.dist = "euclidean", nboot = 10000)
  
  #heatmap of DARs
  #switch row names to standard sample names
  samp.id = sapply(rownames(mat), function(x) {unlist(strsplit(x,":"))[3]})
  inds = match(samp.id, samp.md$sample_id)
  rownames(mat) = samp.md$Sample[inds]
  
  #make column annotation
  pats = sapply(rownames(mat), function(x) {unlist(strsplit(x,"_"))[2]})
  samps = sapply(rownames(mat), function(x) {unlist(strsplit(x,"_"))[4]})
  sample_annotation_df = data.frame(sample = rownames(mat),
                                    grade = sapply(seq(length(rownames(mat))), function(x) {samp.md$Grade[which(samp.md$patient==pats[x] & samp.md$sample_id==samps[x])]}),
                                    ini_rec = sapply(seq(length(rownames(mat))), function(x) {samp.md$ini_rec[which(samp.md$patient==pats[x] & samp.md$sample_id==samps[x])]}),
                                    TMZ = sapply(seq(length(rownames(mat))), function(x) {samp.md$TMZ[which(samp.md$patient==pats[x] & samp.md$sample_id==samps[x])]}),
                                    PCV = sapply(seq(length(rownames(mat))), function(x) {samp.md$PCV[which(samp.md$patient==pats[x] & samp.md$sample_id==samps[x])]}),
                                    RT = sapply(seq(length(rownames(mat))), function(x) {samp.md$Radiotherapy[which(samp.md$patient==pats[x] & samp.md$sample_id==samps[x])]}))
  rownames(sample_annotation_df) = rownames(mat)
  
  mat_plot <- log1p(t(mat))
  mat_plot <- t(scale(t(mat_plot)))
  
  # Sample annotations
  anno_df <- sample_annotation_df[colnames(mat_plot),rev(c("grade", "ini_rec", "TMZ", "PCV", "RT")),drop = FALSE]
  anno_df$grade <- factor(anno_df$grade,levels = c("2", "3", "4"))
  anno_df$ini_rec <- factor(anno_df$ini_rec,levels = c("ini", "rec1","rec2"))
  anno_df$TMZ <- factor(anno_df$TMZ,levels = c("0", "1"))
  anno_df$PCV <- factor(anno_df$PCV,levels = c("0", "1"))
  anno_df$RT <- factor(anno_df$RT,levels = c("0", "1"))
  
  # Annotation bars
  ha <- HeatmapAnnotation(df = anno_df,col = ann_colors,na_col = "grey90")
  
  # Heatmap colour scale
  lim <- 2.2
  col_fun <- colorRamp2(c(-lim, 0, lim),scico::scico(3, palette = "vik"))
  ht <- Heatmap(mat_plot,
                name = "Scaled log\npeak counts",
                col = col_fun,
                heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2)),
                top_annotation = ha,
                cluster_columns = pv_result$hclust,
                cluster_rows = TRUE,
                clustering_method_rows = "ward.D2",
                show_row_names = FALSE,
                column_title = paste("DAR in", cell.type, "cells for", analysis_name))
  
  pdf(paste0("figures_Libra/", cancer, str.mer,"_heatmap_Libra_DAR_", str.ana, "_",grp.backward, "_", grp.forward, "_",cell.type, ".pdf"),
      width = 8,height = 7)
  draw(ht)
  dev.off()
}



########################
####motif enrichment in DARs
########################

#load ArchR object
oligo.ArchR1 = loadArchRProject(path = paste0("~/stageout/ArchR_output_",cancer,"/project1"), force = FALSE, showLogo = TRUE)
#get two main groups of peaks, more and less accessible
#first list is up in forward group, second down in forward group
peaks = list()
peaks[[1]] = results_import$Libra %>%
  filter(cell_type == cell.type, p_val_adj < 0.05, abs(avg_logFC) > 0.5, avg_logFC < 0) %>% pull(gene)
peaks[[2]] = results_import$Libra %>%
  filter(cell_type == cell.type, p_val_adj < 0.05, abs(avg_logFC) > 0.5, avg_logFC > 0) %>% pull(gene)


#do motif enrichment using JASPAR2024 motif set
library(GenomicRanges)
peakSet = lapply(peaks, function(x) {GRanges(x)})

matches <- getMatches(oligo.ArchR1)
motifEnrich <- lapply(peakSet, function(x) {
  customEnrichment(ranges = x, matches = matches)
})


#plot enrichment versus -log10 p val
Nplot = 40
lapply(seq(length(motifEnrich)), function(x) {
  #filter p-adj<0.01
  motifEnrich[[x]] = motifEnrich[[x]][which(motifEnrich[[x]]$mlog10Padj>(-log10(0.01))),]
  #order by enrichment and plot top 40
  df.toplot = motifEnrich[[x]][order(motifEnrich[[x]]$Enrichment, decreasing=T)[1:Nplot],]
  labels = sapply(rownames(df.toplot), function(y) {unlist(strsplit(y,"_"))[1]})
  ggplot(df.toplot, aes(x = Enrichment, y = mlog10Padj, label = labels)) + theme(text = element_text(size = 18)) +
    geom_jitter() + geom_text_repel(size = 4,max.overlaps=20,max.iter = 20000,max.time = 5,force = 5,
                                    force_pull = 0.5, box.padding = 0.5, point.padding = 0.3) + 
    labs(x = "Motif enrichment", y = "-log10(adjusted p value)")
  ggsave(paste0("figures_Libra/",cancer,str.mer,"_motif_enrichment_",str.ana,"_",cell.type,"_",x,".pdf"))
})







