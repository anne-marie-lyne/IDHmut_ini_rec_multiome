library(ArchR)
library(parallel)
library(data.table)
library(readxl)
library(patchwork)
library(viridis)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

#ggplot2 colour palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#set cancer subtype
cancer = "oligo"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)

#minimum number of cells in sample-cell state
#excluded from analysis otherwise
N.min.grp = 100

#ArchR parameters
set.seed(1)
addArchRLocking(locking = T)
addArchRThreads(threads = 1)
addArchRGenome("hg38")

fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_traj_",cancer)
dir.create(fig.dir)

#read in sample list and sample meta data
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))

#load in full ArchR object
oligo.ArchR = loadArchRProject(path = paste0(out.dir,"/project1"), force = T, showLogo = F)


#read in dataframe with trajectory pseudotime
df.traj = fread(paste0("persistent/trajectory/df_file_traj_",cancer,".tsv.gz"), data.table = F)
#add archr style names to df.traj
df.traj$cell.names.archr = sapply(seq(nrow(df.traj)), function(x) {
  cell.n = paste(unlist(strsplit(df.traj$cell[x],"-"))[1:2],collapse="-")
  paste(df.traj$sample[x],cell.n,sep="#")
})

#add sample meta data
df.traj$treatment_PCV = sapply(df.traj$sample, function(x) {samp.md$PCV[which(samp.md$Sample==x)]})
df.traj$treatment_TMZ = sapply(df.traj$sample, function(x) {samp.md$TMZ[which(samp.md$Sample==x)]})
df.traj$treatment_TMZ[df.traj$treatment_TMZ==1] = "TMZ"
df.traj$treatment_TMZ[df.traj$treatment_TMZ==0] = "no_TMZ"
treatments = c("TMZ","PCV","Radiotherapy")
df.traj$treatment_group = sapply(df.traj$sample, function(x) {
  paste(treatments[samp.md[which(samp.md$Sample==x),treatments]==1],collapse="_")
})
df.traj$treatment_group[df.traj$treatment_group==""] = "none"
df.traj$patient = sapply(df.traj$sample, function(x) {unlist(strsplit(x,"_"))[2]})


###############
## TF accessibility
###############


#extract matrix with TF motif deviations from chromVAR
plotVarDev <- getVarDeviations(oligo.ArchR, name = "MotifMatrix", plot = T)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores-astro", width = 5, height = 5, ArchRProj = oligo.ArchR, addDOC = FALSE)

devMat <- getMatrixFromProject(ArchRProj = oligo.ArchR, useMatrix = "MotifMatrix")
zMat <- assay(devMat, "z")

#identify cells in both (fewer cells in ArchR objects due to extra scenicplus filtering and some samples RNA only)
names.both = intersect(colnames(zMat),df.traj$cell.names.archr)
rownames(df.traj) = df.traj$cell.names.archr
df.traj.ss = df.traj[names.both,]
zMat.ss = zMat[,names.both]
#remove all the TFs which are duplicates
zMat.ss = zMat.ss[!duplicated(as.data.frame(as.matrix(zMat.ss))),]


#add column to df.ss with combination of cell state and sample
df.traj.ss$state.samp = paste(df.traj.ss$ICA_score_states,sapply(df.traj.ss$sample, function(x) {unlist(strsplit(x,"_"))[2]}),
                              sapply(df.traj.ss$sample, function(x) {unlist(strsplit(x,"_"))[3]}),sep="_")


######################
## identify TFs specifically active in each state
######################


#anova to find TFs in sample x state, state, sample and patient as explanatory variables
df.to.plot = as.data.frame(as.matrix(t(zMat.ss)))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$sample = df.traj.ss$sample
df.to.plot$patient = df.traj.ss$patient
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","sample","patient"))

#anova and then Tukey test for states
tukey_results.acc <- lapply(tf_cols, function(tf) {
  form <- as.formula(paste0("`", tf, "` ~ ICA_score_states + sample + patient"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$ICA_score_states
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results.acc <- do.call(rbind, tukey_results.acc)


#extract relevant rows
#can do multiple different comparisons...
if(cancer=="oligo"){
  sig.results.opc = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff > 0 & grepl("OPC-like-astro-like",rownames(tukey_results.acc))),]
  sig.results.opc = sig.results.opc[order(sig.results.opc$diff,decreasing=T),]
  sig.results.opc2 = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff > 0 & grepl("OPC-like-intermediate",rownames(tukey_results.acc))),]
  sig.results.opc2 = sig.results.opc2[order(sig.results.opc2$diff,decreasing=T),]
  sig.results.astro = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("OPC-like-astro-like",rownames(tukey_results.acc))),]
  sig.results.astro = sig.results.astro[order(sig.results.astro$diff,decreasing=F),]
  sig.results.astro2 = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("intermediate-astro-like",rownames(tukey_results.acc))),]
  sig.results.astro2 = sig.results.astro2[order(sig.results.astro2$diff,decreasing=F),]
  sig.results.int = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("OPC-like-intermediate",rownames(tukey_results.acc))),]
  sig.results.int = sig.results.int[order(sig.results.int$diff,decreasing=F),]
}else{
  sig.results.opc = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff > 0 & grepl("OPC-like-astro-like",rownames(tukey_results.acc))),]
  sig.results.opc = sig.results.opc[order(sig.results.opc$diff,decreasing=T),]
  sig.results.opc2 = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff > 0 & grepl("OPC-like-intermediate",rownames(tukey_results.acc))),]
  sig.results.opc2 = sig.results.opc2[order(sig.results.opc2$diff,decreasing=T),]
  sig.results.astro = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("OPC-like-astro-like",rownames(tukey_results.acc))),]
  sig.results.astro = sig.results.astro[order(sig.results.astro$diff,decreasing=F),]
  sig.results.astro2 = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("intermediate-astro-like",rownames(tukey_results.acc))),]
  sig.results.astro2 = sig.results.astro2[order(sig.results.astro2$diff,decreasing=F),]
  sig.results.int = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("OPC-like-intermediate",rownames(tukey_results.acc))),]
  sig.results.int = sig.results.int[order(sig.results.int$diff,decreasing=F),]
  sig.results.neuro = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff > 0 & grepl("neuro-like-intermediate",rownames(tukey_results.acc))),]
  sig.results.neuro = sig.results.neuro[order(sig.results.neuro$diff,decreasing=T),]
  sig.results.neuro2 = tukey_results.acc[which(tukey_results.acc$padj < 0.05 & tukey_results.acc$diff < 0 & grepl("OPC-like-neuro-like",rownames(tukey_results.acc))),]
  sig.results.neuro2 = sig.results.neuro2[order(sig.results.neuro2$diff,decreasing=F),]
}


if(cancer=="oligo"){
  TF.rows = unique(c(sig.results.opc$TF[1:8],sig.results.opc2$TF[1:8],sig.results.astro$TF[1:8],sig.results.astro2$TF[1:8]))
  TF.rows = TF.rows[!(TF.rows %in% c("BATF3","BNC2","BATF::JUN","FOSL2::JUN_FOSL2::JUN_MA1130.2","FOSL1::JUN_FOSL1::JUN_MA1128.2","FOSL1"))]
}else{
  TF.rows = unique(c(sig.results.opc$TF[1:8],sig.results.astro$TF[1:8],sig.results.opc2$TF[1:8],sig.results.astro2$TF[1:8],
                     sig.results.neuro2$TF[1:14]))
  #remove some repeated AP-1 factors for visualisation
  TF.rows = TF.rows[!(TF.rows %in% c("BATF3","BNC2","FOSL2::JUN_FOSL2::JUN_MA1130.2","FOS::JUN_FOS::JUN_MA0099.4","FOSL1::JUN_FOSL1::JUN_MA1128.2",
                                     "BATF::JUN","FOSL1::JUNB","FOSL1","ESX1","TAL1::TCF3"))]
}


################
## plot heatmaps of mean value for each state-sample
################

#extract relevant rows
tf.mat <- zMat.ss[TF.rows, , drop = FALSE]
# Group for each maligant state
grp <- factor(df.traj.ss$state.samp)
# Sum within each state/sample group
sum.mat <- t(rowsum(t(tf.mat), group = grp))
# Number of cells in each group
n.per.grp <- table(grp)
# Mean within each group
mean.mat <- sweep(sum.mat, 2, n.per.grp[colnames(sum.mat)], "/")
#keep only thos with at least 50 cells
mean.mat = mean.mat[,n.per.grp>=N.min.grp]

#annotate with grade, treatment metadata etc.
col_anno_df = data.frame(patient = sapply(colnames(mean.mat), function(x) {unlist(strsplit(x,"_"))[2]}))
col_anno_df$ini_rec = sapply(colnames(mean.mat), function(x) {unlist(strsplit(x,"_"))[3]})
col_anno_df$sample = sapply(seq(length(col_anno_df$ini_rec)), function(x) {samp.md$sample_id[which(samp.md$patient==col_anno_df$patient[x] & samp.md$ini_rec==col_anno_df$ini_rec[x])]})
col_anno_df$grade = sapply(col_anno_df$sample, function(x) {samp.md$Grade[which(samp.md$sample_id==x)]})
col_anno_df$TMZ = sapply(col_anno_df$sample, function(x) {samp.md$TMZ[which(samp.md$sample_id==x)]})
col_anno_df$PCV = sapply(col_anno_df$sample, function(x) {samp.md$PCV[which(samp.md$sample_id==x)]})
col_anno_df$RT = sapply(col_anno_df$sample, function(x) {samp.md$Radiotherapy[which(samp.md$sample_id==x)]})
col_anno_df$state = sapply(rownames(col_anno_df), function(x) {unlist(strsplit(x,"_"))[1]})

state.cols = gg_color_hue(4)
names(state.cols) = c("astro-like","intermediate","neuro-like","OPC-like")
cols.list.col = list(grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue"),
                     ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4"),
                     TMZ = c("1" = "grey30", "0" = "gray90"), PCV = c("1" = "grey30", "0" = "gray90"),
                     RT = c("1" = "grey30", "0" = "gray90"),
                     state = state.cols)

top_anno <- HeatmapAnnotation(df = col_anno_df[,c("grade","ini_rec","TMZ","PCV","RT")], col = cols.list.col)
bot_anno <- HeatmapAnnotation(df = col_anno_df[,c("state"),drop = FALSE], col = cols.list.col)

rownames(mean.mat) = sapply(rownames(mean.mat), function(x) {unlist(strsplit(x,"_"))[1]})
col_fun <- colorRamp2(c(-5, 0, 8), c("blue","white","red"))
ht <- Heatmap(
  mean.mat,
  name = "Mean TF motif\ndeviation Z-score",
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = T,
  top_annotation  = top_anno,
  bottom_annotation  = bot_anno,
  column_dend_height = unit(4, "cm")
)
pdf(paste0(fig.dir,"/heatmap_anova_TF_accessibility_",cancer,".pdf"), width = 15, height = 9, useDingbats = FALSE)
draw(ht)
dev.off()


#separately define cells as AP-1 responsive if zMat.ss["FOSL2",]>10

#then try splitting out high AP-1
df.traj.ss$ICA_mes = df.traj.ss$ICA_score_states
df.traj.ss$ICA_mes[which(zMat.ss["FOSL2",]>10)] = "MES"
grp4 <- factor(df.traj.ss$ICA_mes)
sum.mat4 <- t(rowsum(t(tf.mat), group = grp4))
# Number of cells in each group
n.per.grp4 <- table(grp4)
# Mean within each group
mean.mat4 <- sweep(sum.mat4, 2, n.per.grp4[colnames(sum.mat4)], "/")
col_fun <- colorRamp2(c(-3, 0, 5), c("blue","white","red"))
colnames(mean.mat4)[3] = "AP-1-responsive"
rownames(mean.mat4) = sapply(rownames(mean.mat4), function(x) {unlist(strsplit(x,"_"))[1]})
ht <- Heatmap(
  mean.mat4[,],
  name = "Mean TF motif\ndeviation Z-score",
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = T,show_heatmap_legend = T
)
pdf(paste0(fig.dir,"/heatmap_anova_TF_accessibility_ICA_mes_",cancer,".pdf"), width = 4, height = 6, useDingbats = FALSE)
draw(ht)
dev.off()


#compute overlap of AP-1-responsive cells based on ATAC with cells defined as reactive/MES-like based on RNA ICA projection
#(not shown in paper)
#read in mes like cells based on rna
cells.mes.rna = fread(paste0("/mnt/beegfs/home/alyne/stageout/figures_cohort_summary/mes_cells_",cancer,".tsv.gz"), data.table=F)[,1]
cells.mes.atac = df.traj.ss$cell[which(df.traj.ss$ICA_mes=="MES")]
intersect.rna.atac = intersect(cells.mes.rna,cells.mes.atac)
df.traj.ss$mes_rna_atac = sapply(df.traj.ss$cell, function(x) {
  if(x %in% intersect.rna.atac){
    "both"
  }else if(x %in% cells.mes.rna){
    "rna_only"
  }else if(x %in% cells.mes.atac){
    "atac_only"
  }else{
    "other_state"
  }
})
plot.df <- df.traj.ss %>% count(sample, mes_rna_atac) %>% group_by(sample) %>% mutate(prop = n / sum(n))
plot.df$mes_rna_atac = factor(plot.df$mes_rna_atac, levels = c("other_state","rna_only","atac_only","both"))
ggplot(plot.df, aes(x = sample, y = prop, fill = mes_rna_atac)) +
  geom_col(width = 0.8) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample", y = "Proportion of cells", fill = "") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("/mnt/beegfs/home/alyne/stageout/figures_cohort_summary/stacked_barplot_overlap_mes_rna_atac_",cancer,".pdf"))


######################
#if astrocytoma, plot TF accessibility along trajectory for the selected TFs
######################

genes.to.fit.ps1 = unique(c(sig.results.opc$TF[1:8],sig.results.astro$TF[1:8],sig.results.opc2$TF[1:8],sig.results.astro2$TF[1:8]))
genes.to.fit.ps1 = genes.to.fit.ps1[!(genes.to.fit.ps1 %in% c("BATF3","BNC2","FOSL2::JUN_FOSL2::JUN_MA1130.2","FOS::JUN_FOS::JUN_MA0099.4",
                                                              "FOSL1::JUN_FOSL1::JUN_MA1128.2","BATF::JUN","FOSL1::JUNB","FOSL1"))]
genes.to.fit.ps2 = sig.results.neuro2$TF[1:14][which(!(sig.results.neuro2$TF[1:14] %in% c("ESX1","TAL1::TCF3")))]

#mgcv to fit generalised additive model of TF accessibility along pseudotime
library(mgcv)
#branch1
i=1
pseudotime_mat <- matrix(df.traj.ss[,which(colnames(df.traj.ss)==paste0("slingPseudotime_",i))], ncol = 1)
indsNA = which(!is.na(pseudotime_mat))
cellWeights <- matrix(1, nrow = nrow(pseudotime_mat), ncol = 1)


fits.ps1 = list()
for(j in 1:length(genes.to.fit.ps1)){
  print(j)
  df <- data.frame(y = zMat.ss[genes.to.fit.ps1[j], indsNA], pt = pseudotime_mat[indsNA, 1],
                   patient = df.traj.ss$patient[indsNA], sample = df.traj.ss$sample[indsNA])
  df$patient <- factor(df$patient)
  df$sample <- factor(df$sample)
  #fit GAM, with patient and sample as covariates
  fits.ps1[[j]] <- gam(y ~ s(pt, k = 4) + s(patient, bs = "re") + s(sample, bs = "re"),
                       data = df, method = "REML")
}
names(fits.ps1) = genes.to.fit.ps1
saveRDS(fits.ps1, paste0(fig.dir,"/mgcv_fits_ps1_",cancer,".rds"))
#fits.ps1 = readRDS(paste0(fig.dir,"/mgcv_fits_ps1_",cancer,".rds"))


#get 100 equally spaced time points to plot heatmap
pt.grid <- seq(min(pseudotime_mat[indsNA, 1]), max(pseudotime_mat[indsNA, 1]), length.out = 100)
#get average values
heat_mat <- sapply(fits.ps1, function(fit) {
  nd <- data.frame(pt = pt.grid, patient = unique(df.traj.ss$patient)[1], #NB can be anything since going to exclude
                   sample = unique(df.traj.ss$sample)[1])
  predict(fit, newdata = nd, exclude = c("s(patient)", "s(sample)"))
})

heat_mat <- t(heat_mat)
rownames(heat_mat) <- names(fits.ps1)
colnames(heat_mat) <- round(pt.grid, 2)

rownames(heat_mat) = sapply(rownames(heat_mat), function(x) {unlist(strsplit(x,"_"))[1]})
ht = Heatmap(
  heat_mat,
  name = "Smoothed TF motif\ndeviation Z-score",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  # col = col_fun,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(
      x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA)
    )
  }
)
pdf(paste0(fig.dir,"/ps1_heatmap_TFaccessibility_anova_",cancer,".pdf"), width = 10, height = 3, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


#branch2
i=2
pseudotime_mat <- matrix(df.traj.ss[,which(colnames(df.traj.ss)==paste0("slingPseudotime_",i))], ncol = 1)
indsNA = which(!is.na(pseudotime_mat))
cellWeights <- matrix(1, nrow = nrow(pseudotime_mat), ncol = 1)

#fit gams for TFs selcted by anova
fits.ps2 = list()
for(j in 1:length(genes.to.fit.ps2)){
  print(j)
  df <- data.frame(y = zMat.ss[genes.to.fit.ps2[j], indsNA], pt = pseudotime_mat[indsNA, 1],
                   patient = df.traj.ss$patient[indsNA], sample = df.traj.ss$sample[indsNA])
  df$patient <- factor(df$patient)
  df$sample <- factor(df$sample)
  fits.ps2[[j]] <- gam(y ~ s(pt, k = 4) + s(patient, bs = "re") + s(sample, bs = "re"),
                       data = df, method = "REML")
}
names(fits.ps2) = genes.to.fit.ps2
saveRDS(fits.ps2, paste0(fig.dir,"/mgcv_fits_ps2_",cancer,".rds"))
#fits.ps2 = readRDS(paste0(fig.dir,"/mgcv_fits_ps2_",cancer,".rds"))


#get 100 equally spaced time points for plot
pt.grid <- seq(min(pseudotime_mat[indsNA, 1]), max(pseudotime_mat[indsNA, 1]), length.out = 100)
#get average values
heat_mat <- sapply(fits.ps2, function(fit) {
  nd <- data.frame(pt = pt.grid, patient = unique(df.traj.ss$patient)[1], #NB can be anything since going to exclude
                   sample = unique(df.traj.ss$sample)[1])
  predict(fit, newdata = nd, exclude = c("s(patient)", "s(sample)"))
})

heat_mat <- t(heat_mat)
rownames(heat_mat) <- names(fits.ps2)
colnames(heat_mat) <- round(pt.grid, 2)

rownames(heat_mat) = sapply(rownames(heat_mat), function(x) {unlist(strsplit(x,"_"))[1]})
heat_mat = heat_mat[rownames(heat_mat)!="HAND2",]
ht = Heatmap(
  heat_mat,
  name = "Smoothed TF motif\ndeviation Z-score",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  # col = col_fun,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(
      x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA)
    )
  }
)
pdf(paste0(fig.dir,"/ps2_heatmap_TFaccessibility_anova_",cancer,".pdf"), width = 10, height = 2, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()




###############
##TF activity
###############

#more or less same approach, but this time for TF activity
#use gene AUC score as output by SCENIC+ for each TF as activity

out.dir = paste0("/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_",cancer,"/figures_scenicplus_newdb_",cancer)

use.extended = F
#read in region AUC scores
TF.AUCg.d = fread(paste0(out.dir,"/direct_gene_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
if(use.extended){
  TF.AUCg.e = fread(paste0(out.dir,"/extended_gene_based_AUC_",cancer,".tsv.gz"), header=T, data.table = F)
}else{
  TF.AUCg.e = data.frame(Cell = TF.AUCg.d$Cell)
}
TF.AUCg = cbind(TF.AUCg.d,TF.AUCg.e[,-1])
rownames(TF.AUCg) = TF.AUCg$Cell
rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {unlist(strsplit(x,"___",fixed=T))[1]})
TF.AUCg$Cell = NULL
rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {
  parts = unlist(strsplit(x,"-",fixed=T))
  paste(parts[3],paste0(parts[1],"-",parts[2]),sep="#")
})
TF.AUCg.ss = TF.AUCg[names.both,]

#only keep +/+
TF.AUCg.ss = TF.AUCg.ss[,grepl("+/+",colnames(TF.AUCg.ss),fixed=T)]


#anova to find TFs in sample x state, for activity
df.to.plot = as.data.frame(as.matrix(TF.AUCg.ss))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$sample = df.traj.ss$sample
df.to.plot$patient = df.traj.ss$patient
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","sample","patient"))

tukey_results.act <- lapply(tf_cols, function(tf) {
  form <- as.formula(paste0("`", tf, "` ~ ICA_score_states + sample + patient"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$ICA_score_states
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results.act <- do.call(rbind, tukey_results.act)


#extract relevant rows
if(cancer=="oligo"){
  sig.results.opc = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff > 0 & grepl("OPC-like-astro-like",rownames(tukey_results.act))),]
  sig.results.opc = sig.results.opc[order(sig.results.opc$diff,decreasing=T),]
  sig.results.opc2 = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff > 0 & grepl("OPC-like-intermediate",rownames(tukey_results.act))),]
  sig.results.opc2 = sig.results.opc2[order(sig.results.opc2$diff,decreasing=T),]
  sig.results.astro = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("OPC-like-astro-like",rownames(tukey_results.act))),]
  sig.results.astro = sig.results.astro[order(sig.results.astro$diff,decreasing=F),]
  sig.results.astro2 = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("intermediate-astro-like",rownames(tukey_results.act))),]
  sig.results.astro2 = sig.results.astro2[order(sig.results.astro2$diff,decreasing=F),]
  sig.results.int = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("OPC-like-intermediate",rownames(tukey_results.act))),]
  sig.results.int = sig.results.int[order(sig.results.int$diff,decreasing=F),]
}else{
  sig.results.opc = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff > 0 & grepl("OPC-like-astro-like",rownames(tukey_results.act))),]
  sig.results.opc = sig.results.opc[order(sig.results.opc$diff,decreasing=T),]
  sig.results.opc2 = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff > 0 & grepl("OPC-like-intermediate",rownames(tukey_results.act))),]
  sig.results.opc2 = sig.results.opc2[order(sig.results.opc2$diff,decreasing=T),]
  sig.results.astro = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("OPC-like-astro-like",rownames(tukey_results.act))),]
  sig.results.astro = sig.results.astro[order(sig.results.astro$diff,decreasing=F),]
  sig.results.astro2 = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("intermediate-astro-like",rownames(tukey_results.act))),]
  sig.results.astro2 = sig.results.astro2[order(sig.results.astro2$diff,decreasing=F),]
  sig.results.int = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("OPC-like-intermediate",rownames(tukey_results.act))),]
  sig.results.int = sig.results.int[order(sig.results.int$diff,decreasing=F),]
  sig.results.neuro = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff > 0 & grepl("neuro-like-intermediate",rownames(tukey_results.act))),]
  sig.results.neuro = sig.results.neuro[order(sig.results.neuro$diff,decreasing=T),]
  sig.results.neuro2 = tukey_results.act[which(tukey_results.act$padj < 0.05 & tukey_results.act$diff < 0 & grepl("OPC-like-neuro-like",rownames(tukey_results.act))),]
  sig.results.neuro2 = sig.results.neuro2[order(sig.results.neuro2$diff,decreasing=F),]
}


#add some known literature TFs, or TFs which were specifically accessible to visualise
if(cancer=="oligo"){
  TF.cols = unique(c(sig.results.opc$TF[1:10],sig.results.astro$TF[1:8],sig.results.opc2$TF[1:10],sig.results.astro2$TF[1:8],
                     "FOSB_direct_+/+_(68g)","NFIA_direct_+/+_(654g)",
                     "FOS_direct_+/+_(106g)","JUNB_direct_+/+_(75g)","JUN_direct_+/+_(149g)"))
}else{
  TF.cols = unique(c(sig.results.opc$TF[1:6],sig.results.astro$TF[1:6],sig.results.neuro$TF[1:8],sig.results.opc2$TF[1:6],sig.results.astro2$TF[1:6],
                     "FOSB_direct_+/+_(71g)","FOSL2_direct_+/+_(265g)","FOS_direct_+/+_(119g)", "ID4_direct_+/+_(278g)",sig.results.neuro2$TF[1:2],
                     "JUNB_direct_+/+_(83g)","JUN_direct_+/+_(116g)","NFIA_direct_+/+_(1143g)"))
}

#extract relevant columns
tf.mat <- t(TF.AUCg.ss[, TF.cols, drop = FALSE])
# Scale each TF across cells (accessiblity was already z-scores)
tf.mat.scaled <- t(scale(t(tf.mat)))
# Group for each malignant state
grp <- factor(df.traj.ss$state.samp)
# Mean scaled value within each state/sample group
sum.mat <- t(rowsum(t(tf.mat.scaled), group = grp))
n.per.grp <- table(grp)
mean.mat.scaled <- sweep(sum.mat, 2, n.per.grp[colnames(sum.mat)], "/")
#keep only thos with at least 100 cells
mean.mat.scaled <- mean.mat.scaled[, n.per.grp >= N.min.grp]
rownames(mean.mat.scaled) = sapply(rownames(mean.mat.scaled), function(x) {unlist(strsplit(x,"_"))[1]})

#heatmap annotation, treatment, grade etc
col_anno_df = data.frame(patient = sapply(colnames(mean.mat.scaled), function(x) {unlist(strsplit(x,"_"))[2]}))
col_anno_df$ini_rec = sapply(colnames(mean.mat.scaled), function(x) {unlist(strsplit(x,"_"))[3]})
col_anno_df$sample = sapply(seq(length(col_anno_df$ini_rec)), function(x) {samp.md$sample_id[which(samp.md$patient==col_anno_df$patient[x] & samp.md$ini_rec==col_anno_df$ini_rec[x])]})
col_anno_df$grade = sapply(col_anno_df$sample, function(x) {samp.md$Grade[which(samp.md$sample_id==x)]})
col_anno_df$TMZ = sapply(col_anno_df$sample, function(x) {samp.md$TMZ[which(samp.md$sample_id==x)]})
col_anno_df$PCV = sapply(col_anno_df$sample, function(x) {samp.md$PCV[which(samp.md$sample_id==x)]})
col_anno_df$RT = sapply(col_anno_df$sample, function(x) {samp.md$Radiotherapy[which(samp.md$sample_id==x)]})
col_anno_df$state = sapply(rownames(col_anno_df), function(x) {unlist(strsplit(x,"_"))[1]})

state.cols = gg_color_hue(4)
names(state.cols) = c("astro-like","intermediate","neuro-like","OPC-like")
cols.list.col = list(grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue"),
                     ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4"),
                     TMZ = c("1" = "grey30", "0" = "gray90"),
                     PCV = c("1" = "grey30", "0" = "gray90"),
                     RT = c("1" = "grey30", "0" = "gray90"),
                     state = state.cols)

top_anno <- HeatmapAnnotation(df = col_anno_df[,c("grade","ini_rec","TMZ","PCV","RT")], col = cols.list.col)
bot_anno <- HeatmapAnnotation(df = col_anno_df[,c("state"),drop = FALSE], col = cols.list.col)

ht <- Heatmap(
  mean.mat.scaled,
  name = "Mean scaled\nTF activity",
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = T,
  top_annotation  = top_anno,
  bottom_annotation  = bot_anno,
  column_dend_height = unit(4, "cm")
)
pdf(paste0(fig.dir,"/heatmap_anova_TF_activity_scaled_",cancer,".pdf"), width = 15, height = 9, useDingbats = FALSE)
draw(ht)
dev.off()


#can plot histogram to see distribution of AP-1 factor used for AP1-responsiveness categorisation
# ggplot(as.data.frame(as.matrix(TF.AUCg.ss)), aes(x = .data[["JUN_direct_+/+_(149g)"]])) + geom_histogram(bins=100)
# ggsave(paste0(fig.dir,"/hist_JUN_activity_",cancer,".pdf"),width=10)



#splitting out high AP-1 based on JUN activity since this is in both scenicplus outputs
df.traj.ss$ICA_mes = df.traj.ss$ICA_score_states
if(cancer=="astro"){
  df.traj.ss$ICA_mes[which(TF.AUCg.ss[,"JUN_direct_+/+_(116g)"]>0.04)] = "MES"
}else{
  df.traj.ss$ICA_mes[which(TF.AUCg.ss[,"JUN_direct_+/+_(149g)"]>0.035)] = "MES"
}

grp4 <- factor(df.traj.ss$ICA_mes)
tf.mat.scaled <- t(scale(t(tf.mat)))
sum.mat4 <- t(rowsum(t(tf.mat.scaled), group = grp4))
# Number of cells in each group
n.per.grp4 <- table(grp4)
# Mean within each group
mean.mat4 <- sweep(sum.mat4, 2, n.per.grp4[colnames(sum.mat4)], "/")
col_fun <- colorRamp2(c(-1, 0, 2), c("blue","white","red"))
colnames(mean.mat4)[3] = "AP-1-responsive"
rownames(mean.mat4) = sapply(rownames(mean.mat4), function(x) {unlist(strsplit(x,"_"))[1]})
ht <- Heatmap(
  mean.mat4[,],
  name = "Mean scaled\nTF activity",
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = T,show_heatmap_legend = T
)
pdf(paste0(fig.dir,"/heatmap_anova_TF_activity_scaled_ICA_mes_",cancer,".pdf"), width = 3.5, height = 7, useDingbats = FALSE)
draw(ht)
dev.off()




#if astrocytoma, plot along trajectory
#fit GAMs as above for accessibility

genes.to.fit.ps1 = unique(c(sig.results.opc$TF[1:6],sig.results.astro$TF[1:6],sig.results.opc2$TF[1:6],
                            sig.results.astro2$TF[1:6],"FOSB_direct_+/+_(71g)","FOSL2_direct_+/+_(265g)","FOS_direct_+/+_(119g)","ID4_direct_+/+_(278g)",
                            "JUNB_direct_+/+_(83g)","JUN_direct_+/+_(116g)","NFIA_direct_+/+_(1143g)"))
genes.to.fit.ps2 = unique(c(sig.results.neuro$TF[1:8],sig.results.neuro2$TF[1:2]))

library(mgcv)
i=1
#branch1
pseudotime_mat <- matrix(df.traj.ss[,which(colnames(df.traj.ss)==paste0("slingPseudotime_",i))], ncol = 1)
indsNA = which(!is.na(pseudotime_mat))
cellWeights <- matrix(1, nrow = nrow(pseudotime_mat), ncol = 1)

fits.ps1 = list()
for(j in 1:length(genes.to.fit.ps1)){
  print(j)
  df <- data.frame(y = TF.AUCg.ss[indsNA,genes.to.fit.ps1[j]], pt = pseudotime_mat[indsNA, 1],
                   patient = df.traj.ss$patient[indsNA], sample = df.traj.ss$sample[indsNA])
  df$patient <- factor(df$patient)
  df$sample <- factor(df$sample)
  
  fits.ps1[[j]] <- gam(y ~ s(pt, k = 4) + s(patient, bs = "re") + s(sample, bs = "re"),
                       data = df, method = "REML")
}
names(fits.ps1) = genes.to.fit.ps1
saveRDS(fits.ps1, paste0(fig.dir,"/mgcv_act_fits_ps1_",cancer,".rds"))
#fits.ps1 = readRDS(paste0(fig.dir,"/mgcv_act_fits_ps1_",cancer,".rds"))


#get 100 equally spaced time points
pt.grid <- seq(min(pseudotime_mat[indsNA, 1]), max(pseudotime_mat[indsNA, 1]), length.out = 100)
#get average values
heat_mat <- sapply(fits.ps1, function(fit) {
  nd <- data.frame(pt = pt.grid, patient = unique(df.traj.ss$patient)[1], #NB can be anything since going to exclude
                   sample = unique(df.traj.ss$sample)[1])
  predict(fit, newdata = nd, exclude = c("s(patient)", "s(sample)"))
})

heat_mat <- t(heat_mat)
rownames(heat_mat) <- names(fits.ps1)
colnames(heat_mat) <- round(pt.grid, 2)
#scale
heat_mat <- t(scale(t(heat_mat)))

rownames(heat_mat) = sapply(rownames(heat_mat), function(x) {unlist(strsplit(x,"_"))[1]})
ht = Heatmap(
  heat_mat,
  name = "Scaled smoothed\nTF activity",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  # col = col_fun,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(
      x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA)
    )
  }
)
pdf(paste0(fig.dir,"/ps1_heatmap_TFactivity_anova_",cancer,".pdf"), width = 10, height = 3, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#branch2
i=2
pseudotime_mat <- matrix(df.traj.ss[,which(colnames(df.traj.ss)==paste0("slingPseudotime_",i))], ncol = 1)
indsNA = which(!is.na(pseudotime_mat))
cellWeights <- matrix(1, nrow = nrow(pseudotime_mat), ncol = 1)

#fit gams for TFs selcted by anova
fits.ps2 = list()
for(j in 1:length(genes.to.fit.ps2)){
  print(j)
  df <- data.frame(y = TF.AUCg.ss[indsNA,genes.to.fit.ps2[j]], pt = pseudotime_mat[indsNA, 1],
                   patient = df.traj.ss$patient[indsNA], sample = df.traj.ss$sample[indsNA])
  df$patient <- factor(df$patient)
  df$sample <- factor(df$sample)
  fits.ps2[[j]] <- gam(y ~ s(pt, k = 4) + s(patient, bs = "re") + s(sample, bs = "re"),
                       data = df, method = "REML")
}
names(fits.ps2) = genes.to.fit.ps2
saveRDS(fits.ps2, paste0(fig.dir,"/mgcv_act_fits_ps2_",cancer,".rds"))
#fits.ps2 = readRDS(paste0(fig.dir,"/mgcv_act_fits_ps2_",cancer,".rds"))


#get 100 equally spaced time points for plot
pt.grid <- seq(min(pseudotime_mat[indsNA, 1]), max(pseudotime_mat[indsNA, 1]), length.out = 100)
#get average values
heat_mat <- sapply(fits.ps2, function(fit) {
  nd <- data.frame(pt = pt.grid,patient = unique(df.traj.ss$patient)[1], #NB can be anything since going to exclude
                   sample = unique(df.traj.ss$sample)[1])
  predict(fit, newdata = nd, exclude = c("s(patient)", "s(sample)"))
})

heat_mat <- t(heat_mat)
rownames(heat_mat) <- names(fits.ps2)
colnames(heat_mat) <- round(pt.grid, 2)
#scale
heat_mat <- t(scale(t(heat_mat)))

rownames(heat_mat) = sapply(rownames(heat_mat), function(x) {unlist(strsplit(x,"_"))[1]})
heat_mat = heat_mat[rownames(heat_mat)!="TFDP2",]
heat_mat = heat_mat[rownames(heat_mat)!="HOXA9",]
ht = Heatmap(
  heat_mat,
  name = "Scaled smoothed\nTF activity",
  cluster_rows = T,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  # col = col_fun,
  rect_gp = gpar(col = NA),
  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(
      x, y,
      w + unit(0.3, "pt"),
      h + unit(0.3, "pt"),
      gp = gpar(fill = fill, col = NA)
    )
  }
)
pdf(paste0(fig.dir,"/ps2_heatmap_TFactivity_anova_",cancer,".pdf"), width = 10, height = 2, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()




