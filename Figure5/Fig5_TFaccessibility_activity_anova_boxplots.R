library(ArchR)
library(parallel)
library(data.table)
library(readxl)
library(patchwork)
library(viridis)
library(dplyr)
library(tidyr)
library(ggpubr)


#code to produce boxplots in Figure 5
#based on anaova analysis of accessibility and activity by grade and treatment status in IDH-A


#code to identify DARs, make heatmaps and do motif enrichment is in folder: Libra_ATAC_DAR_analysis



#define subtype
cancer = "astro"
base.dir = "/mnt/beegfs/home/alyne"
setwd(base.dir)
out.dir = paste0("stageout/ArchR_output_",cancer)
dir.create(out.dir)

set.seed(1)

addArchRLocking(locking = T)
addArchRThreads(threads = 1)
addArchRGenome("hg38")

fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_traj_",cancer)
dir.create(fig.dir)

#read in sample list
samps = fread(paste0("stagein/scenicplus_",cancer,"/sample_list.txt"),data.table = F, header = F)[,1]

oligo.ArchR = loadArchRProject(path = paste0(out.dir,"/project1"), force = T, showLogo = F)

#read in dataframe with trajectory pseudotime value for each cell
df.traj = fread(paste0("persistent/trajectory/df_file_traj_",cancer,".tsv.gz"), data.table = F)
#add archr style names to df.traj
df.traj$cell.names.archr = sapply(seq(nrow(df.traj)), function(x) {
  cell.n = paste(unlist(strsplit(df.traj$cell[x],"-"))[1:2],collapse="-")
  paste(df.traj$sample[x],cell.n,sep="#")
})


#read in sample meta data
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))
df.traj$treatment_PCV = sapply(df.traj$sample, function(x) {samp.md$PCV[which(samp.md$Sample==x)]})
df.traj$treatment_TMZ = sapply(df.traj$sample, function(x) {samp.md$TMZ[which(samp.md$Sample==x)]})
df.traj$treatment_TMZ[df.traj$treatment_TMZ==1] = "TMZ"
df.traj$treatment_TMZ[df.traj$treatment_TMZ==0] = "no_TMZ"
treatments = c("TMZ","PCV","Radiotherapy")
df.traj$treatment_group = sapply(df.traj$sample, function(x) {
  paste(treatments[samp.md[which(samp.md$Sample==x),treatments]==1],collapse="_")
})
df.traj$treatment_group[df.traj$treatment_group==""] = "none"


#extract chromVAR deviations
plotVarDev <- getVarDeviations(oligo.ArchR, name = "MotifMatrix", plot = T)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores-astro", width = 5, height = 5, ArchRProj = oligo.ArchR, addDOC = FALSE)

devMat <- getMatrixFromProject(ArchRProj = oligo.ArchR, useMatrix = "MotifMatrix")
zMat <- assay(devMat, "z")

#extract cells which are in both (fwer cells in chromVAR due to scenic+ filtering and three samples without ATAC)
names.both = intersect(colnames(zMat),df.traj$cell.names.archr)
rownames(df.traj) = df.traj$cell.names.archr
df.traj.ss = df.traj[names.both,]
zMat.ss = zMat[,names.both]


######
##Motif accessibility
######


#first look at impact of grade using only untreated samples in IDH-A
#use intermediate cells as otherwise have problem with very variables numbers of cells in each sample/grade
#grade and sample as covariates to control for sample-specific behaviour
df.to.plot = as.data.frame(as.matrix(t(zMat.ss)))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$grade = factor(df.traj.ss$grade,levels = c(2,3,4))
df.to.plot$treatment = df.traj.ss$treatment
df.to.plot$sample = df.traj.ss$sample
df.to.plot = df.to.plot[which(df.to.plot$ICA_score_states=="intermediate"),]
df.to.plot = df.to.plot[which(df.to.plot$treatment=="untreated"),]
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","grade","sample","treatment"))

tukey_results <- lapply(tf_cols, function(tf) {
  
  form <- as.formula(paste0("`", tf, "` ~ grade + sample"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$grade
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results <- do.call(rbind, tukey_results)

sig.results.pos = subset(tukey_results, padj < 0.05 & diff > 0)
sig.results.neg = subset(tukey_results, padj < 0.05 & diff < 0)

head(sig.results.pos[order(sig.results.pos$diff,decreasing=T),],50)
head(sig.results.neg[order(sig.results.neg$diff,decreasing=F),],50)



#plot boxplot of accessibility by grade

cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "blue3")
gr.genes = c("POU3F2","GSC","OTX2","FEZF2","NFIA","NFIX_NFIX_MA0671.2","DBP","RFX5")

#version showing all cells
treatment.genes = c(gr.genes)
df.list = list()
for(i in 1:length(treatment.genes)){
  df3 = data.frame(accessibility = zMat.ss[treatment.genes[i],], state = df.traj.ss[["ICA_score_states"]], sample = df.traj.ss[["sample"]],
                   treatment = df.traj.ss[["treatment"]], grade = df.traj.ss[["grade"]])
  df3 = df3[which(df3$treatment=="untreated"),]
  df3 = df3[which(df3$state=="intermediate"),]
  df3$TF = treatment.genes[i]
  df.list[[i]] = df3
}
df = do.call("rbind",df.list)
df$grade = factor(df$grade, levels = c(2,3,4))
df$TF[which(df$TF=="NFIX_NFIX_MA0671.2")] = "NFIX_NFIX"
df$TF = factor(df$TF, levels = c("POU3F2","GSC","OTX2","FEZF2","NFIA","NFIX_NFIX","DBP","RFX5"))
ggplot(df, aes(x=grade, y=accessibility, fill=grade)) + facet_wrap(~TF, ncol = length(treatment.genes)) +
  stat_compare_means(aes(group = grade),comparisons = comparisons,method = "wilcox.test",
                     label = "p.signif", step.increase = 0.1, hide.ns = F) + scale_fill_manual(values = cols.grade) +
  geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),size = 0.005,alpha = 0.5) +
  labs(title="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+ ylim(-10, 17)
ggsave(paste0(fig.dir,"/boxplot_grade_untreated_all_TFs_accessibility_up_down_",cancer,".pdf"),width=7,height=5)







#next treatment
#anova of accessiblility TMZ-treated compared to untreated in grade 4 IDH-A
df.to.plot = as.data.frame(as.matrix(t(zMat.ss)))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$grade = df.traj.ss$grade
df.to.plot$treatment = df.traj.ss$treatment_group
df.to.plot$treatment[which(df.to.plot$treatment=="TMZ_Radiotherapy")] = "TMZ"
df.to.plot$sample = df.traj.ss$sample
df.to.plot = df.to.plot[which(df.to.plot$ICA_score_states=="intermediate"),]
df.to.plot = df.to.plot[which(df.to.plot$grade==4),]
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","grade","sample","treatment"))

tukey_results <- lapply(tf_cols, function(tf) {
  
  form <- as.formula(paste0("`", tf, "` ~ treatment + sample"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$treatment
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results <- do.call(rbind, tukey_results)


sig.results.pos = subset(tukey_results, padj < 0.05 & diff > 0)
sig.results.neg = subset(tukey_results, padj < 0.05 & diff < 0)

head(sig.results.pos[order(sig.results.pos$diff,decreasing=T),],50)
head(sig.results.neg[order(sig.results.neg$diff,decreasing=F),],50)



#boxplot

treatment.genes = c("EGR1","FOSL2","FOS::JUND","BACH1","PITX1","GSC","OTX2","FEZF2")
df.list = list()
for(i in 1:length(treatment.genes)){
  df3 = data.frame(accessibility = zMat.ss[treatment.genes[i],], state = df.traj.ss[["ICA_score_states"]], sample = df.traj.ss[["sample"]],
                   treatment = df.traj.ss[["treatment_TMZ"]], grade = df.traj.ss[["grade"]])
  df3 = df3[which(df3$grade==4),]
  df3 = df3[which(df3$state=="intermediate"),]
  df3$TF = treatment.genes[i]
  df.list[[i]] = df3
}
df = do.call("rbind",df.list)
df$treatment = factor(df$treatment, levels = c("no_TMZ","TMZ"))
df$TF = factor(df$TF, levels = treatment.genes)
ggplot(df, aes(x=TF, y=accessibility, fill=treatment)) + 
  geom_boxplot(position = position_dodge(width = 0.8),alpha = 0.5,outlier.shape = NA) + labs(title="") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),size = 0.005,alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(-12, 50) +
  stat_compare_means(aes(group = treatment),method = "wilcox.test",label = "p.signif")
ggsave(paste0(fig.dir,"/boxplot_all_treatment_gr4_TFs_accessibility_",cancer,".pdf"),width=7,height=5)




########
##TF activity
########


out.dir = paste0("/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_",cancer,"/figures_scenicplus_newdb_",cancer)

#read in SCENIC+ gene AUC score for each motif in each cell
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
TF.AUCg.ss = t(TF.AUCg.ss[,grepl("+/+",colnames(TF.AUCg.ss),fixed=T)])




#do simlar analysis to above but with activity



#TMZ-treated grade 4 IDH-A compared to untreated
df.to.plot = as.data.frame(as.matrix(t(TF.AUCg.ss)))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$grade = df.traj.ss$grade
df.to.plot$treatment = df.traj.ss$treatment_group
df.to.plot$treatment[which(df.to.plot$treatment=="TMZ_Radiotherapy")] = "TMZ"
df.to.plot$sample = df.traj.ss$sample
df.to.plot = df.to.plot[which(df.to.plot$ICA_score_states=="intermediate"),]
df.to.plot = df.to.plot[which(df.to.plot$grade==4),]
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","grade","sample","treatment"))

tukey_results <- lapply(tf_cols, function(tf) {
  
  form <- as.formula(paste0("`", tf, "` ~ treatment + sample"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$treatment
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results <- do.call(rbind, tukey_results)

sig.results.pos = subset(tukey_results, padj < 0.05 & diff > 0)
sig.results.neg = subset(tukey_results, padj < 0.05 & diff < 0)

head(sig.results.pos[order(sig.results.pos$diff,decreasing=T),],50)
head(sig.results.neg[order(sig.results.neg$diff,decreasing=F),],50)



#plot boxplot of activity in intermediate cells
treatment.genes = c("EGR1_direct_+/+_(119g)","FOSL2_direct_+/+_(265g)","JUN_direct_+/+_(116g)",
                    "POU3F3_direct_+/+_(25g)","SOX4_direct_+/+_(260g)","MEIS2_direct_+/+_(43g)","ASCL1_direct_+/+_(171g)",
                    "CUX2_direct_+/+_(200g)","SOX8_direct_+/+_(169g)")
df.list = list()
for(i in 1:length(treatment.genes)){
  df3 = data.frame(accessibility = TF.AUCg.ss[treatment.genes[i],], state = df.traj.ss[["ICA_score_states"]], sample = df.traj.ss[["sample"]],
                   treatment = df.traj.ss[["treatment_TMZ"]], grade = df.traj.ss[["grade"]])
  df3 = df3[which(df3$grade==4),]
  df3 = df3[which(df3$state=="intermediate"),]
  df3$TF = treatment.genes[i]
  df.list[[i]] = df3
}
df = do.call("rbind",df.list)
df$treatment = factor(df$treatment, levels = c("no_TMZ","TMZ"))
df$TF = factor(df$TF, levels = treatment.genes)
ggplot(df, aes(x=TF, y=accessibility, fill=treatment)) +
  stat_compare_means(aes(group = treatment),method = "wilcox.test",label = "p.signif") + 
  geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),size = 0.005,alpha = 0.5) +
  labs(title="",y = "Activity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(fig.dir,"/boxplot_all_treatment_gr4_TFs_activity_",cancer,".pdf"),width=7,height=5)



#for early response genes and AP-1, plot activity in astro-like cells as we see something different here
gr.genes = c("EGR1_direct_+/+_(119g)","FOSL2_direct_+/+_(265g)","JUN_direct_+/+_(116g)")

df.list = list()
for(i in 1:length(gr.genes)){
  df3 = data.frame(accessibility = TF.AUCg.ss[gr.genes[i],], state = df.traj.ss[["ICA_score_states"]], sample = df.traj.ss[["sample"]],
                   treatment = df.traj.ss[["treatment_TMZ"]], grade = df.traj.ss[["grade"]])
  df3 = df3[which(df3$grade==4),]
  df3 = df3[which(df3$state=="astro-like"),]
  df3$TF = gr.genes[i]
  df.list[[i]] = df3
}
df = do.call("rbind",df.list)
df$treatment = factor(df$treatment, levels = c("no_TMZ","TMZ"))
df$TF = factor(df$TF, levels = gr.genes)
ggplot(df, aes(x=TF, y=accessibility, fill=treatment)) +
  stat_compare_means(aes(group = treatment),method = "wilcox.test",label = "p.signif") + 
  geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),size = 0.005,alpha = 0.5) +
  labs(title="",y = "Activity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(fig.dir,"/boxplot_all_treatment_gr4_TFs_activity_",cancer,"_astrolike.pdf"),width=3,height=5)





#next impact of grade in untreated IDH-A samples
#TF activity, grade and sample as covariates
df.to.plot = as.data.frame(as.matrix(t(TF.AUCg.ss)))
df.to.plot$ICA_score_states = df.traj.ss$ICA_score_states
df.to.plot$grade = factor(df.traj.ss$grade,levels = c(2,3,4))
df.to.plot$treatment = df.traj.ss$treatment
df.to.plot$sample = df.traj.ss$sample
df.to.plot = df.to.plot[which(df.to.plot$ICA_score_states=="intermediate"),]
df.to.plot = df.to.plot[which(df.to.plot$treatment=="untreated"),]
tf_cols <- setdiff(colnames(df.to.plot), c("ICA_score_states","grade","sample","treatment"))

tukey_results <- lapply(tf_cols, function(tf) {
  
  form <- as.formula(paste0("`", tf, "` ~ grade + sample"))
  
  fit <- aov(form, data = df.to.plot)
  
  tk <- TukeyHSD(fit)$grade
  
  data.frame(TF = tf,comparison = rownames(tk),diff = tk[, "diff"],padj = tk[, "p adj"])
})
tukey_results <- do.call(rbind, tukey_results)

sig.results.pos = subset(tukey_results, padj < 0.05 & diff > 0)
sig.results.neg = subset(tukey_results, padj < 0.05 & diff < 0)

head(sig.results.pos[order(sig.results.pos$diff,decreasing=T),],50)
head(sig.results.neg[order(sig.results.neg$diff,decreasing=F),],50)




cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "blue3")
gr.genes = c("POU3F3_direct_+/+_(25g)","SOX11_direct_+/+_(75g)","SOX4_direct_+/+_(260g)",
             "VAX2_direct_+/+_(16g)","ID4_direct_+/+_(278g)","NFIA_direct_+/+_(1143g)",
             "NR3C2_direct_+/+_(240g)","RFX4_direct_+/+_(544g)")

comparisons <- list(c("2", "3"),c("2", "4"),c("3", "4"))
df.list = list()
for(i in 1:length(gr.genes)){
  df3 = data.frame(accessibility = TF.AUCg.ss[gr.genes[i],], state = df.traj.ss[["ICA_score_states"]], sample = df.traj.ss[["sample"]],
                   treatment = df.traj.ss[["treatment"]], grade = df.traj.ss[["grade"]])
  df3 = df3[which(df3$treatment=="untreated"),]
  df3 = df3[which(df3$state=="intermediate"),]
  df3$TF = gr.genes[i]
  df.list[[i]] = df3
}
df1 = do.call("rbind",df.list)
df1$grade = factor(df1$grade, levels = c(2,3,4))
df1$TF = factor(df1$TF, levels = gr.genes)
ggplot(df1, aes(x=grade, y=accessibility, fill = grade)) + facet_wrap(~TF, ncol = length(gr.genes))+
  stat_compare_means(aes(group = grade),comparisons = comparisons,method = "wilcox.test",
    label = "p.signif", step.increase = 0.2, hide.ns = F)  + lims(y = c(-0.02,0.48)) +
  geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA) + scale_fill_manual(values = cols.grade) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),size = 0.005,alpha = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.8),alpha = 0.5,outlier.shape = NA) + labs(title="",y = "Activity") +
  theme(strip.text = element_text(size = 6))
ggsave(paste0(fig.dir,"/boxplot_all_grade_untreated_TFs_activity_",cancer,".pdf"),width=10,height=4)




#finally, score cells for the relevant ICA metacomponents (stress/AP-1, Reactive/MES-like and neuro-like)
#plot boxplot of scores for IDH-A grade 4 grouped by whether they were TMZ-treated or not
#facetted by maligannt cell state

data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
fig.str = "malig"
if(cancer=="astro"){
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
}else{
  oligo.ss = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
}
oligo.ss.filt = subset(oligo.ss, subset = grade==4)



treatment.genes = c("ICA_3_1","ICA_8_1","ICA_13_1")
IC.names = c("Stress/AP-1","Reactive/MES-like","Neuro-like")
plot.list = list()
for(i in 1:length(treatment.genes)){
  df3 = data.frame(score = oligo.ss.filt@meta.data[,treatment.genes[i]], state = oligo.ss.filt@meta.data[,"ICA_score_states"], sample = oligo.ss.filt@meta.data[,"orig.ident"],
                   treatment = oligo.ss.filt@meta.data[,"treatment"], grade = oligo.ss.filt@meta.data[,"grade"])
  df3$TF = treatment.genes[i]
  df3$state = factor(df3$state)
  df3$treatment[which(df3$treatment=="untreated")] = "no_TMZ"
  df3$treatment[which(df3$treatment=="treated")] = "TMZ"
  df3$treatment = factor(df3$treatment,levels = c("no_TMZ","TMZ"))
  plot.list[[i]] = ggplot(df3, aes(treatment, score, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.3, size = 0.2) +
  facet_wrap(~state,ncol=4) + labs(title=IC.names[i],y = "Score",x = "") +
  stat_compare_means(method = "wilcox.test",comparisons = list(c("no_TMZ", "TMZ")),label = "p.signif") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x=element_blank())
  if(i%%3!=0){
    plot.list[[i]] = plot.list[[i]] + guides(fill="none")
  }
}
wrap_plots(plot.list)
ggsave(paste0(fig.dir,"/boxplot_all_treatment_gr4_ICA_scores_",cancer,"_by_state.pdf"),width=12,height=4)


