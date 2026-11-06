library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Signac)
library(readxl)
library(tidyr)
library(ggpubr)
library(rstatix)
library(textshape)
library(compositions)
library(easyCODA)
library(vegan)
library(data.table)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
dir.create(fig.dir)


#try again getting rid of endo, myeloid, opc, lymphoid


#load in Seurat objects for each subtype, both all cells and only malignant
astro = readRDS("scenicplus_astro/oligo.ss_integrated_RPCA_all_astro_doub_rem.rds")
oligo = readRDS("scenicplus_oligo/oligo.ss_integrated_RPCA_all_oligo_doub_rem.rds")
astro.malig = readRDS("scenicplus_astro/oligo.ss_integrated_RPCA_malig_astro_doub_rem.rds")
oligo.malig = readRDS("scenicplus_oligo/oligo.ss_integrated_RPCA_malig_oligo.rds")


#rename mac to myeloid
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Macrophage")] = "Myeloid"


#sICA thresholds so can read in relevant metagene info for each subtype

cancer = "astro"
cor.thresh = 0.45
fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.astro = colnames(df.annot)

cancer = "oligo"
cor.thresh = 0.5
fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.oligo = colnames(df.annot)



#find top projected dimension for each cell
oligo.malig$top_ICA_oligo = sapply(seq(length(colnames(oligo.malig))), function(x) {prog.orig.order.oligo[which.max(oligo.malig@assays$ICA@layers$data[,x])]})
astro.malig$top_ICA_astro = sapply(seq(length(colnames(astro.malig))), function(x) {prog.orig.order.astro[which.max(astro.malig@assays$ICA@layers$data[,x])]})


#redefine malignant cells as reactive/MES if their top dimension is Reactive astrocyte
oligo.malig$ICA_mes = oligo.malig$ICA_score_states
oligo.malig$ICA_mes[which(oligo.malig$top_ICA_oligo=="Reactive/MES-like")] = "reactive"

astro.malig$ICA_mes = astro.malig$ICA_score_states
astro.malig$ICA_mes[which(astro.malig$top_ICA_astro=="Reactive/MES-like")] = "reactive"


#read in myeloid cell states
# fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
# anno = fread(paste0(fig.dir,"/cell_annotations_minus_S0009.tsv.gz"), data.table=F)
fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_myeloid'
anno = fread(paste0(fig.dir,"/cell_annotations_all_samples.tsv.gz"), data.table=F)
colnames(anno) = c("cell.name","orig.ident", "status","patient", "grade", "sample", "disease", "cluster" ,"cluster_annotation")

#add patient and sample id
astro$patient = sapply(colnames(astro), function(x) {unlist(strsplit(x,"_"))[2]})
oligo$patient = sapply(colnames(oligo), function(x) {unlist(strsplit(x,"_"))[2]})
astro$sample = sapply(colnames(astro), function(x) {unlist(strsplit(x,"_"))[4]})
oligo$sample = sapply(colnames(oligo), function(x) {unlist(strsplit(x,"_"))[4]})


#combine endo and VLMC
astro$cell.type.sr.cl.imm.postmerge[which(astro$iCNV.malig.ident=="malignant")] = paste0("Malignant_",astro$ICA_score_states[which(astro$iCNV.malig.ident=="malignant")])
names.mes = names(astro.malig$ICA_mes[which(astro.malig$ICA_mes=="reactive")])
inds.mes = match(names.mes, names(astro$orig.ident))
astro$cell.type.sr.cl.imm.postmerge[inds.mes] = "Malignant_reactive"

#replace Macrophage with lineage-based annotation
inds.astro = match(anno$cell.name, names(astro$patient))
astro$cell.type.sr.cl.imm.postmerge[inds.astro[which(!is.na(inds.astro))]] = paste0("Myeloid_",anno$cluster_annotation[which(!is.na(inds.astro))])
#replace cells which are cycling with annotation
cells.cyc = names(astro.malig$orig.ident)[which(astro.malig$cyc.status=="cycling")]

inds.astro = match(cells.cyc, names(astro$patient))
astro$cell.type.sr.cl.imm.postmerge[inds.astro] = "Malignant_cycling"

astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Neuron"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "Neuron"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Glial"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Astro")] = "Glial"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="B-cells")] = "Lymphoid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="T-cells")] = "Lymphoid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="NK cells")] = "Lymphoid"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endo"
astro$cell.type.sr.cl.imm.postmerge[which(astro$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"


oligo$cell.type.sr.cl.imm.postmerge[which(oligo$iCNV.malig.ident=="malignant")] = paste0("Malignant_",oligo$ICA_score_states[which(oligo$iCNV.malig.ident=="malignant")])
names.mes = names(oligo.malig$ICA_mes[which(oligo.malig$ICA_mes=="reactive")])
inds.mes = match(names.mes, names(oligo$orig.ident))
oligo$cell.type.sr.cl.imm.postmerge[inds.mes] = "Malignant_reactive"
#replace Macrophage with lineage-based annotation
inds.oligo = match(anno$cell.name, names(oligo$patient))
oligo$cell.type.sr.cl.imm.postmerge[inds.oligo[which(!is.na(inds.oligo))]] = paste0("Myeloid_",anno$cluster_annotation[which(!is.na(inds.oligo))])
#replace cells which are cycling with annotation
cells.cyc = names(oligo.malig$orig.ident)[which(oligo.malig$cyc.status=="cycling")]
inds.oligo = match(cells.cyc, names(oligo$patient))
oligo$cell.type.sr.cl.imm.postmerge[inds.oligo] = "Malignant_cycling"

oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Glutamatergic")] = "Neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="GABAergic")] = "Neuron"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Oligo")] = "Glial"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Astro")] = "Glial"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="B-cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="T-cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="NK cells")] = "Lymphoid"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="VLMC")] = "Endo"
oligo$cell.type.sr.cl.imm.postmerge[which(oligo$cell.type.sr.cl.imm.postmerge=="Endo")] = "Endothelial"


fwrite(as.data.frame(oligo$cell.type.sr.cl.imm.postmerge),
       paste0(fig.dir,"/cell_types_oligo_postmerge_separate.tsv.gz"), row.names = T)
fwrite(as.data.frame(astro$cell.type.sr.cl.imm.postmerge),
       paste0(fig.dir,"/cell_types_astro_postmerge_separate.tsv.gz"), row.names = T)


fwrite(as.data.frame(oligo$cell.type.sr.cl.imm.postmerge),
       paste0(fig.dir,"/cell_types_oligo_postmerge.tsv.gz"), row.names = T)
fwrite(as.data.frame(astro$cell.type.sr.cl.imm.postmerge),
       paste0(fig.dir,"/cell_types_astro_postmerge.tsv.gz"), row.names = T)
cell.annot.oligo = fread(paste0(fig.dir,"/cell_types_oligo_postmerge.tsv.gz"),data.table=F, header=T)
cell.annot.astro = fread(paste0(fig.dir,"/cell_types_astro_postmerge.tsv.gz"),data.table=F, header=T)
# cell.annot.oligo[which(cell.annot.oligo[,2]=="Myeloid_Macrophage SPP1+"),2] = "Myeloid_Macrophage"
# cell.annot.astro[which(cell.annot.astro[,2]=="Myeloid_Macrophage SPP1+"),2] = "Myeloid_Macrophage"
oligo$cell.type.sr.cl.imm.postmerge = cell.annot.oligo[,2]
astro$cell.type.sr.cl.imm.postmerge = cell.annot.astro[,2]


df.to.print.oligo = oligo@meta.data[,c("orig.ident","sample","grade","patient","ini_rec",
                                       "ICA_score_states","cell.type.sr.cl.imm.postmerge","iCNV.malig.ident")]
colnames(df.to.print.oligo) = c("sample_id","sample","grade","patient","ini_rec",
                                "malignant_lineage_states","full_cell_states","malignant")
df.to.print.astro = astro@meta.data[,c("orig.ident","sample","grade","patient","ini_rec",
                                       "ICA_score_states","cell.type.sr.cl.imm.postmerge","iCNV.malig.ident")]
colnames(df.to.print.astro) = c("sample_id","sample","grade","patient","ini_rec",
                                "malignant_lineage_states","full_cell_states","malignant")
fwrite(df.to.print.oligo,
       paste0(fig.dir,"/cell_metadata_oligo.tsv.gz"), row.names = T)
fwrite(df.to.print.astro,
       paste0(fig.dir,"/cell_metadata_astro.tsv.gz"), row.names = T)


#NB shouldn't be any of the filtered cell types other than some myeloid cells which were QC filtered
df.astro.pat = data.frame(cell.type = astro$cell.type.sr.cl.imm.postmerge, ini_rec = astro$ini_rec,
                          patient = astro$patient, cancer = "Astrocytoma", sample = astro$sample)
df.astro.pat <- df.astro.pat %>%
  filter(!cell.type %in% c("Lymphoid","Endothelial","Myeloid","OPC"))
df_counts.astro <- df.astro.pat %>%
  count(patient, ini_rec, cancer, cell.type, sample, name = "count") %>%
  group_by(sample) %>%
  complete(cell.type = unique(df.astro.pat$cell.type),
           fill = list(count = 0)) %>%
  fill(patient, ini_rec, cancer, .direction = "downup") %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()


df.oligo.pat = data.frame(cell.type = oligo$cell.type.sr.cl.imm.postmerge, ini_rec = oligo$ini_rec,
                          patient = oligo$patient, cancer = "Oligodendroglioma", sample = oligo$sample)
df.oligo.pat <- df.oligo.pat %>%
  filter(!cell.type %in% c("Lymphoid","Endothelial","Myeloid","OPC"))
df_counts.oligo <- df.oligo.pat %>%
  count(patient, ini_rec, cancer, cell.type, sample, name = "count") %>%
  group_by(sample) %>%
  complete(cell.type = unique(df.oligo.pat$cell.type),
           fill = list(count = 0)) %>%
  fill(patient, ini_rec, cancer, .direction = "downup") %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

df_counts = rbind(df_counts.astro,df_counts.oligo)
df_counts$cancer = factor(df_counts$cancer)


#transformation for statistical tests
gli.counts0 <- df_counts %>%
  select(sample, cell.type, count) %>%
  pivot_wider(names_from = cell.type, values_from = count, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()
gli.counts = gli.counts0 +1

### --------------------------------------------------------------------------------
### chiPower function
chiPower<- function(data, close=TRUE, power=1, chi=TRUE, BoxCox=TRUE, center=TRUE) {
  # function to compute chiPower transformation with various options
  # close: close data
  # power: tranformation (i.e., value of lambda)
  # chi: chi-square standardization
  # BoxCox: apply 1/power rescaling (without the subtraction of 1)
  # center: center the final result by column means
  foo <- data^power
  if(close)   foo <- foo / rowSums(foo)
  if(chi)     foo <- sweep(foo, 2, sqrt(colMeans(foo)), FUN="/")
  if(BoxCox)  foo <- (1/power) * foo
  if(center)  foo <- sweep(foo, 2, colMeans(foo))
  if(!center) foo <- foo - 1/sqrt(ncol(data))
  data.chiPower <- foo
  data.chiPower
}
### ----------------------------------------------------------------------------------

### logratio distances (with column weight constant in this case, i.e. weight=FALSE)
### Use CLRs to get logratio distances
gli.clr <- CLR(CLOSE(gli.counts), weight=FALSE)$LR
### LRA is PCA of CLRs (or you can get it directly using function LRA)
gli.lra <- PCA(gli.clr, weight=FALSE)
### Principal row coordinates define logratio geometry (last dimension has no variance)
gli.lra.rpc <- gli.lra$rowpcoord[,-(ncol(gli.counts)-1)]
### Logratio distances from CLRs
gli.clr.dist <- dist(gli.clr)/sqrt(ncol(gli.counts))


### loop on powers from 0.001 to 1 and find optimal using Procrustes
### comparison with chiPower transforms in two different ways
### 1. Spearman rank correlation between distances
### 2. Procrustes correlation between geometries
npower <- 100
gli.cor <- rep(0,npower)
gli.pro <- rep(0,npower)  
### Original composition on data set with 0s
gli0.comp <- CLOSE(gli.counts0)    
### start of loop (takes a while) -----------
for(power in 1:npower) {
  foo.chiPower <- chiPower(gli0.comp, power=power/npower)
  foo.dist     <- dist(foo.chiPower)
  foo.rpc      <- PCA(foo.chiPower, weight=FALSE)$rowpcoord[,1:(ncol(gli0.comp)-1)]
  gli.cor[power] <- cor(gli.clr.dist, foo.dist, method="spearman")
  gli.pro[power] <- protest(gli.lra.rpc, foo.rpc, permutations=0)$t0
}


fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_cohort_summary'
### plot the Procrustes (black) and Spearman correlations (blue)
pdf(paste0(fig.dir,"/chiPower_correlations_all_malig_combined_full_clusters.pdf"))
plot(gli.pro, type="l", lwd=2, xaxt="n", bty="n", ylim=c(0.3,1),
     xlab="Power of chiPower transformation", ylab="Procrustes correlation / Spearman")
axis(1, at=seq(0,100,10), labels=seq(0, 1, 0.1))
### blue curve for Spearman correlation
lines(gli.cor, col="blue")
### (Procrustes correlation maximum of 0.901 at power=0.25
###  Spearman correlation maximum of 0.753 at power=0.30)
dev.off()


#transform counts using Procrustes max
gli.chiPower <- chiPower(gli0.comp, power=max(gli.pro))
fwrite(gli.chiPower, paste0(fig.dir,"/chiPower_all_samples_merge_cell_types.tsv.gz"))




#get subtype-specific 
gli.chiPower.astro = gli.chiPower[rownames(gli.chiPower) %in% unique(astro$sample),]
gli.chiPower.oligo = gli.chiPower[rownames(gli.chiPower) %in% unique(oligo$sample),]


#add transformed values to df_counts
df_counts <- df_counts %>%
  mutate(chiPower = gli.chiPower[cbind(sample, cell.type)])


#read in meta.data
samp.md = as.data.frame(read_excel("~/persistent/sample_summary_astro_oligo_anon.xlsx"))
samp.md$DOB <- as.Date(samp.md$DOB, origin = "1899-12-30")
samp.md$Surgery_date <- as.Date(samp.md$Surgery_date, origin = "1899-12-30")
samp.md$DDN <- as.Date(samp.md$DDN, origin = "1899-12-30")


#use only initial samples and predict survival
#start with astro since all patients have died


##########
## IDH-A
##########

#find relevant samples to extract
row.inds = match(unique(astro$sample), samp.md$sample_id)
samples.ini = unique(astro$sample)[which(samp.md$ini_rec[row.inds]=="ini")]

#extract survival time for these patients
row.inds.ini = match(samples.ini, samp.md$sample_id)
patients = samp.md$patient[row.inds.ini]

surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")])/365
})

gli.chiPower.astro.ini = gli.chiPower.astro[samples.ini,]


#fit linear models for each cell type/state
lm.list = lapply(colnames(gli.chiPower.astro.ini), function(x) {
  df <- data.frame(
    surv.time = surv.time,
    predictor = gli.chiPower.astro.ini[, x]
  )
  lm(surv.time ~ predictor, data = df)
})
names(lm.list) <- colnames(gli.chiPower.astro.ini)

lapply(lm.list, function(x) {summary(x)})
#malignant cycling only thing significant


######
##same thing from recurrence samples
samples.rec1 = unique(astro$sample)[which(samp.md$ini_rec[row.inds]=="rec1")]

#extract survival time for these patients
row.inds.rec1 = match(samples.rec1, samp.md$sample_id)
patients = samp.md$patient[row.inds.rec1]

surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="rec1")])/365
})

gli.chiPower.astro.rec1 = gli.chiPower.astro[samples.rec1,]


#fit linear models for each cell type/state
lm.list = lapply(colnames(gli.chiPower.astro.rec1), function(x) {
  df <- data.frame(
    surv.time = surv.time,
    predictor = gli.chiPower.astro.rec1[, x]
  )
  lm(surv.time ~ predictor, data = df)
})
names(lm.list) <- colnames(gli.chiPower.astro.rec1)

lapply(lm.list, function(x) {summary(x)})

#astro-like significant, cycling close to being significant




#######
### IDH-O
#######

library(survival)
library(survminer)

#for oligo, have censored data, so need cox model
#find relevant samples to extract
row.inds = match(unique(oligo$sample), samp.md$sample_id)
samples.ini = unique(oligo$sample)[which(samp.md$ini_rec[row.inds]=="ini")]

#extract survival time for these patients
row.inds.ini = match(samples.ini, samp.md$sample_id)
patients = samp.md$patient[row.inds.ini]

surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")])/365
})

status = sapply(patients, function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==x)]
  state = rel.rows[!is.na(rel.rows)]
  if(state=="Dead"){
    1
  }else{
    0
  }
})

gli.chiPower.oligo.ini = gli.chiPower.oligo[samples.ini,]

cox.list <- lapply(colnames(gli.chiPower.oligo.ini), function(x) {
  df <- data.frame(
    surv.time = surv.time,
    status = status,
    predictor = gli.chiPower.oligo.ini[, x]
  )
  coxph(Surv(surv.time, status) ~ predictor, data = df)
  
})
names(cox.list) <- colnames(gli.chiPower.oligo.ini)

cox.results <- do.call(rbind, lapply(names(cox.list), function(x) {
  s <- summary(cox.list[[x]])
  data.frame(celltype = x,HR = exp(coef(s)),coef = coef(s),p.value = s$coefficients[,"Pr(>|z|)"],concordance = s$concordance[1])
}))
cox.results[order(cox.results$p.value), ]


#cycling and reactive significant



#same but from recurrence samples
row.inds = match(unique(oligo$sample), samp.md$sample_id)
samples.rec1 = unique(oligo$sample)[which(samp.md$ini_rec[row.inds]=="rec1")]

#extract survival time for these patients
row.inds.rec1 = match(samples.rec1, samp.md$sample_id)
patients = samp.md$patient[row.inds.rec1]

surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="rec1")])/365
})

status = sapply(patients, function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==x)]
  state = rel.rows[!is.na(rel.rows)]
  if(state=="Dead"){
    1
  }else{
    0
  }
})

gli.chiPower.oligo.rec1 = gli.chiPower.oligo[samples.rec1,]

#one model per cell type
cox.list <- lapply(colnames(gli.chiPower.oligo.rec1), function(x) {
  df <- data.frame(
    surv.time = surv.time,
    status = status,
    predictor = gli.chiPower.oligo.rec1[, x]
  )
  coxph(Surv(surv.time, status) ~ predictor, data = df)
})
names(cox.list) <- colnames(gli.chiPower.oligo.rec1)

cox.results <- do.call(rbind, lapply(names(cox.list), function(x) {
  s <- summary(cox.list[[x]])
  data.frame(celltype = x,HR = exp(coef(s)),coef = coef(s),p.value = s$coefficients[,"Pr(>|z|)"],concordance = s$concordance[1])
}))
cox.results[order(cox.results$p.value), ]

#cycling the only thing significant





#scatter plots, Figure 7a

row.inds = match(unique(astro$sample), samp.md$sample_id)
samples.ini = unique(astro$sample)
#extract survival time for these patients
row.inds.ini = match(samples.ini, samp.md$sample_id)
patients = samp.md$patient[row.inds.ini]
ini_rec = samp.md$ini_rec[row.inds.ini]
surv.time = sapply(seq(length(patients)), function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==patients[x])]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==patients[x] & samp.md$ini_rec==ini_rec[x])])/365
})
gli.chiPower.astro.ini = gli.chiPower.astro[samples.ini,]

df = data.frame(survival_time = surv.time, ini_rec = ini_rec,
                cycling = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_cycling"),"prop"]})),
                astro.like = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_astro-like"),"prop"]})),
                reactive = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_reactive"),"prop"]})),
                inflammatory = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Inflammatory"),"prop"]})),
                macrophage = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Macrophage"),"prop"]})),
                neuron = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Neuron"),"prop"]}))
)
cell.types <- setdiff(colnames(df), c("survival_time","ini_rec"))
plots <- lapply(cell.types, function(cell) {
  ggplot(df, aes_string(x = cell, y = "survival_time")) + geom_point(size = 3, aes(color=ini_rec)) +
    geom_smooth(method = "lm", se = TRUE) +labs(x = cell,y = "Survival time (years)") +theme_classic() + 
    scale_color_manual(values = c("ini" = "violet", "rec1" = "purple")) +theme(text = element_text(size = 16))
})
wrap_plots(plots, ncol = 3)
ggsave(paste0(fig.dir,"/scatter_astro_ini_and_rec_survival.pdf"),width = 15, height = 6)


lapply(cell.types, function(cell) {
  cor.test(df[df$ini_rec=="ini","survival_time"],df[df$ini_rec=="ini",cell])
})
lapply(cell.types, function(cell) {
  cor.test(df[df$ini_rec=="rec1","survival_time"],df[df$ini_rec=="rec1",cell])
})






library(survival)
library(survminer)

row.inds = match(unique(oligo$sample), samp.md$sample_id)
samples.ini = unique(oligo$sample)[which(samp.md$ini_rec[row.inds]=="ini")]
#extract survival time for these patients
row.inds.ini = match(samples.ini, samp.md$sample_id)
patients = samp.md$patient[row.inds.ini]
surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="ini")])/365
})
status = sapply(patients, function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==x)]
  state = rel.rows[!is.na(rel.rows)]
  if(state=="Dead"){
    1
  }else{
    0
  }
})
gli.chiPower.oligo.ini = gli.chiPower.oligo[samples.ini,]
#make dataframe of proportions and survival time
df = data.frame(survival_time = surv.time, 
                status = status,
                cycling = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_cycling"),"prop"]})),
                astro.like = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_astro-like"),"prop"]})),
                reactive = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_reactive"),"prop"]})),
                inflammatory = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Inflammatory"),"prop"]})),
                macrophage = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Macrophage"),"prop"]})),
                neuron = unlist(sapply(samples.ini, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Neuron"),"prop"]}))
)
cell.types <- setdiff(colnames(df), "survival_time")
df$cycling_group <- ifelse(df$cycling > median(df$cycling),"High cycling","Low cycling")
df$cycling_group <- factor(df$cycling_group,levels = c("Low cycling", "High cycling"),labels = c("Low", "High"))
fit <- survfit(Surv(survival_time, status) ~ cycling_group,data = df)
fit$strata <- setNames(fit$strata,sub("^cycling_group=", "", names(fit$strata)))
p <- surv_pvalue(fit, data = df)$pval
pdf(paste0(fig.dir,"/km_oligo_ini_survival_cycling.pdf"),width = 6,height = 5)
plt = ggsurvplot(fit, data = df, pval = paste0("Log-rank test p = ", signif(p, 2)),risk.table = TRUE,
           legend.title = "group",title = "Cycling",break.time.by = 5,xlim = c(0, 25))
plt$plot <- plt$plot +
  theme(text = element_text(size = 16), plot.title = element_text(size = 18), axis.title = element_text(size = 16),
    axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
print(plt)
dev.off()


df$reactive_group <- ifelse(df$reactive > median(df$reactive),"High reactive","Low reactive")
df$status = status
df$reactive_group <- factor(df$reactive_group,levels = c("Low reactive", "High reactive"),labels = c("Low", "High"))
fit <- survfit(Surv(survival_time, status) ~ reactive_group,data = df)
fit$strata <- setNames(fit$strata,sub("^reactive_group=", "", names(fit$strata)))
p <- surv_pvalue(fit, data = df)$pval
pdf(paste0(fig.dir,"/km_oligo_ini_survival_reactive.pdf"),width = 6,height = 5)
plt = ggsurvplot(fit, data = df, pval = paste0("Log-rank test p = ", signif(p, 2)),risk.table = TRUE,
                 legend.title = "group",title = "Reactive/MES-like",break.time.by = 5,xlim = c(0, 25))
plt$plot <- plt$plot +
  theme(text = element_text(size = 16), plot.title = element_text(size = 18), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
print(plt)
dev.off()


row.inds = match(unique(oligo$sample), samp.md$sample_id)
samples.rec1 = unique(oligo$sample)[which(samp.md$ini_rec[row.inds]=="rec1")]
#extract survival time for these patients
row.inds.rec1 = match(samples.rec1, samp.md$sample_id)
patients = samp.md$patient[row.inds.rec1]
surv.time = sapply(patients, function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==x)]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==x & samp.md$ini_rec=="rec1")])/365
})
status = sapply(patients, function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==x)]
  state = rel.rows[!is.na(rel.rows)]
  if(state=="Dead"){
    1
  }else{
    0
  }
})
gli.chiPower.oligo.rec1 = gli.chiPower.oligo[samples.rec1,]
#make dataframe of proportions and survival time
df = data.frame(survival_time = surv.time, 
                cycling = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_cycling"),"prop"]})),
                astro.like = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_astro-like"),"prop"]})),
                reactive = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Malignant_reactive"),"prop"]})),
                inflammatory = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Inflammatory"),"prop"]})),
                macrophage = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Myeloid_Macrophage"),"prop"]})),
                neuron = unlist(sapply(samples.rec1, function(x) {df_counts[which(df_counts$sample==x & df_counts$cell.type=="Neuron"),"prop"]}))
)
cell.types <- setdiff(colnames(df), "survival_time")


df$cycling_group <- ifelse(df$cycling > median(df$cycling),"High cycling","Low cycling")
df$status = status
df$cycling_group <- factor(df$cycling_group,levels = c("Low cycling", "High cycling"),labels = c("Low", "High"))
fit <- survfit(Surv(survival_time, status) ~ cycling_group,data = df)
fit$strata <- setNames(fit$strata,sub("^cycling_group=", "", names(fit$strata)))
p <- surv_pvalue(fit, data = df)$pval
pdf(paste0(fig.dir,"/km_oligo_rec1_survival_cycling.pdf"),width = 6,height = 5)
plt = ggsurvplot(fit, data = df, pval = paste0("Log-rank test p = ", signif(p, 2)),risk.table = TRUE,
                 legend.title = "group",title = "Cycling",break.time.by = 5,xlim = c(0, 25))
plt$plot <- plt$plot +
  theme(text = element_text(size = 16), plot.title = element_text(size = 18), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
print(plt)
dev.off()


df$reactive_group <- ifelse(df$reactive > median(df$reactive),"High reactive","Low reactive")
df$status = status
df$reactive_group <- factor(df$reactive_group,levels = c("Low reactive", "High reactive"),labels = c("Low", "High"))
fit <- survfit(Surv(survival_time, status) ~ reactive_group,data = df)
fit$strata <- setNames(fit$strata,sub("^reactive_group=", "", names(fit$strata)))
p <- surv_pvalue(fit, data = df)$pval
pdf(paste0(fig.dir,"/km_oligo_rec1_survival_reactive.pdf"),width = 6,height = 5)
plt = ggsurvplot(fit, data = df, pval = paste0("Log-rank test p = ", signif(p, 2)),risk.table = TRUE,
                 legend.title = "group",title = "Reactive/MES-like",break.time.by = 5,xlim = c(0, 25))
plt$plot <- plt$plot +
  theme(text = element_text(size = 16), plot.title = element_text(size = 18), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
print(plt)
dev.off()






#make correlation plots of chiPower values within cancer

library(ComplexHeatmap)
library(circlize)
M = cor(gli.chiPower.astro, method="spearman")
col_fun <- colorRamp2(c(-0.8, 0, 0.8), c("blue","white","red"))
ht = Heatmap(as.matrix(M),
             name = "Correlation",
             col = col_fun, 
             # row_split = 2, column_split = 2,
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2",
             row_names_max_width = unit(8.5, "cm"),
             row_names_gp = gpar(fontsize = 16),
             column_names_gp = gpar(fontsize = 16))
pdf(paste0(fig.dir,"/correlation_spearman_chiPower_astro_full_clusters_merged.pdf"), width = 12, height = 11, useDingbats = FALSE)
draw(ht,heatmap_legend_side = "right",padding = unit(c(30, 10, 10, 20), "mm"))
dev.off()


gli.chiPower.oligo = gli.chiPower.oligo[,which(colnames(gli.chiPower.oligo)!="Malignant_neuro-like")]
M = cor(gli.chiPower.oligo, method="spearman")
M[is.na(M)] = 0
col_fun <- colorRamp2(c(-0.8, 0, 0.8), c("blue","white","red"))
ht = Heatmap(as.matrix(M),
             name = "Correlation",
             col = col_fun, 
             # row_split = 2, column_split = 2,
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2",
             row_names_max_width = unit(8.5, "cm"),
             row_names_gp = gpar(fontsize = 16),
             column_names_gp = gpar(fontsize = 16),
             heatmap_legend_param = list(
               title_gp = gpar(fontsize = 18),
               labels_gp = gpar(fontsize = 16),
               legend_height = unit(5, "cm")
             ),)
pdf(paste0(fig.dir,"/correlation_spearman_chiPower_oligo_full_clusters_merged.pdf"), width = 12, height = 11, useDingbats = FALSE)
draw(ht,heatmap_legend_side = "right",padding = unit(c(30, 10, 10, 20), "mm"))
dev.off()






#pca analysis of cell type/state proportions
to_pca = gli.chiPower[which(rownames(gli.chiPower) %in% samp.md$sample_id[which(samp.md$cancer=="astro")]),]
pca.astro = prcomp(to_pca)

#make data frame with meta data about samples
meta <- data.frame(sample = rownames(gli.chiPower)) %>%
  left_join(samp.md, by = c("sample" = "sample_id")) %>%
  mutate(treatment = case_when(
    TMZ == 1 & PCV == 0 ~ "TMZ",
    PCV == 1 & TMZ == 0 ~ "PCV",
    TRUE ~ "none"))

pca_df.astro <- data.frame(PC1 = pca.astro$x[,1], PC2 = pca.astro$x[,2], sample = rownames(pca.astro$x)) %>%
  left_join(meta, by = "sample")


to.plot = c("ini_rec","Grade","treatment")
plot.list = list()
for(i in 1:3){
  plot.list[[i]] = ggplot(pca_df.astro, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
    scale_color_brewer(palette = "Set1") + theme_classic()
  if(i==1){
    cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
    plot.list[[i]] = ggplot(pca_df.astro, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
      scale_color_manual(values = cols.ini_rec) + theme_classic()
  }
  if(i==2){
    cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
    plot.list[[i]] = ggplot(pca_df.astro, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
      scale_color_manual(values = cols.grade) + theme_classic()
  }
}
wrap_plots(plot.list, ncol=3)
ggsave(paste0(fig.dir,"/pca_compositional_astro_full_clusters.pdf"), width = 12, height = 4)


pca_df.astro$treatment_comb = sapply(pca_df.astro$treatment, function(x) {
  if(x %in% c("PCV","TMZ")){
    "treated"
  }else{
    "none"
  }
})
#what drives first two PCA axes
lm1 <- lm(PC1 ~ treatment_comb + Grade + ini_rec, data = pca_df.astro)
lm2 <- lm(PC2 ~ treatment_comb + Grade + ini_rec, data = pca_df.astro)
summary(lm1)
summary(lm2)

loadings <- as.data.frame(pca.astro$rotation[,1:2])
loadings$celltype <- rownames(loadings)
loadings[order(abs(loadings$PC1), decreasing = TRUE), ]
loadings[order(abs(loadings$PC2), decreasing = TRUE), ]



#fit linear models predicting transformed cell proportions from grade, treatment, ini_rec status
lm_df.astro <- to_pca %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% left_join(meta, by = "sample")
lm_df.astro$treatment_comb = sapply(lm_df.astro$treatment, function(x) {
  if(x %in% c("PCV","TMZ")){
    "treated"
  }else{
    "none"
  }
})

#sort out column names as problem with dashes etc
colnames(lm_df.astro) = gsub("-","",colnames(lm_df.astro))
colnames(lm_df.astro) = gsub(" ","",colnames(lm_df.astro))
colnames(lm_df.astro) = gsub("+","",colnames(lm_df.astro),fixed=T)
lm_df.astro$treatment_comb <- factor(lm_df.astro$treatment_comb, levels = c("none", "treated"))
features <- c("Glial","Malignant_astrolike","Malignant_cycling","Malignant_intermediate","Malignant_neurolike",       
              "Malignant_OPClike","Malignant_reactive","Myeloid_Inflammatory","Myeloid_Macrophage",
              "Myeloid_MacrophageSPP1","Myeloid_Microglia","Neuron")
results <- lapply(features, function(x) {
  fit <- lm(reformulate(c("Grade", "treatment_comb", "ini_rec"), response = x), data = lm_df.astro)
  
  co <- summary(fit)$coefficients
  ci <- confint(fit)
  
  data.frame(feature = x, grade_beta = co["Grade", "Estimate"], grade_p = co["Grade", "Pr(>|t|)"], grade_low = ci["Grade", 1],
    grade_high = ci["Grade", 2], treatment_beta = co["treatment_combtreated", "Estimate"],
    treatment_p = co["treatment_combtreated", "Pr(>|t|)"], treatment_low = ci["treatment_combtreated", 1],
    treatment_high = ci["treatment_combtreated", 2])
})
results <- do.call(rbind, results)

plot_df <- results %>%
  select(feature,grade_beta, grade_p, grade_low, grade_high,
         treatment_beta, treatment_p, treatment_low, treatment_high) %>%
  pivot_longer(cols = -feature, names_to = c("variable", ".value"), names_pattern = "(grade|treatment)_(.*)")
plot_df$feature = gsub("like","-like",plot_df$feature)
plot_df$feature = gsub("SPP1","SPP1+",plot_df$feature)
plot_df$feature = factor(plot_df$feature, 
                         levels = c(unique(plot_df$feature)[!grepl("^M",unique(plot_df$feature))], 
                                    unique(plot_df$feature)[grepl("Myeloid_",unique(plot_df$feature))],
                                    unique(plot_df$feature)[grepl("Malignant_",unique(plot_df$feature))]))
#plot figure of coefficient estimates and confidence intervals
ggplot(plot_df, aes(x = beta, y = feature, color = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar( aes(xmin = low, xmax = high), orientation = "y", height = 0,
    position = position_dodge(width = 0.35)) +
  geom_point( size = 3, position = position_dodge(width = 0.35)) +
  theme_classic() + labs( x = "Regression coefficient", y = NULL, color = NULL)
ggsave(paste0(fig.dir,"/lm_estimates_astro_full_clusters.pdf"), width = 6, height = 5)





###same for oligos
to_pca = gli.chiPower[which(rownames(gli.chiPower) %in% samp.md$sample_id[which(samp.md$cancer=="oligo")]),]
pca.oligo = prcomp(to_pca)

pca_df.oligo <- data.frame(PC1 = pca.oligo$x[,1], PC2 = pca.oligo$x[,2], sample = rownames(pca.oligo$x)) %>%
  left_join(meta, by = "sample")

to.plot = c("ini_rec","Grade","treatment")
plot.list = list()
for(i in 1:3){
  plot.list[[i]] = ggplot(pca_df.oligo, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
    scale_color_brewer(palette = "Set1") + theme_classic()
  if(i==1){
    cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
    plot.list[[i]] = ggplot(pca_df.oligo, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
      scale_color_manual(values = cols.ini_rec) + theme_classic()
  }
  if(i==2){
    cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
    plot.list[[i]] = ggplot(pca_df.oligo, aes(PC1, PC2)) + geom_point(aes(color = factor(.data[[to.plot[i]]]))) +
      scale_color_manual(values = cols.grade) + theme_classic()
  }
}
wrap_plots(plot.list, ncol=3)
ggsave(paste0(fig.dir,"/pca_compositional_oligo_full_clusters.pdf"), width = 12, height = 3)

#what drives first two PCA axes
lm1 <- lm(PC1 ~ treatment + Grade + ini_rec, data = pca_df.oligo)
lm2 <- lm(PC2 ~ treatment + Grade + ini_rec, data = pca_df.oligo)
summary(lm1)
summary(lm2)

loadings <- as.data.frame(pca.oligo$rotation[,1:2])
loadings$celltype <- rownames(loadings)
loadings[order(abs(loadings$PC1), decreasing = TRUE), ]




#fit linear models predicting transformed cell proportions from grade, treatment, ini_rec status
lm_df.oligo <- to_pca %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% left_join(meta, by = "sample")
lm_df.oligo$treatment_comb = sapply(lm_df.oligo$treatment, function(x) {
  if(x %in% c("PCV","TMZ")){
    "treated"
  }else{
    "none"
  }
})

#sort out column names as problem with dashes etc
colnames(lm_df.oligo) = gsub("-","",colnames(lm_df.oligo))
colnames(lm_df.oligo) = gsub(" ","",colnames(lm_df.oligo))
colnames(lm_df.oligo) = gsub("+","",colnames(lm_df.oligo),fixed=T)
lm_df.oligo$treatment <- factor(lm_df.oligo$treatment, levels = c("none", "PCV", "TMZ"))
lm_df.oligo$treatment_comb <- factor(lm_df.oligo$treatment_comb, levels = c("none", "treated"))
features <- c("Glial","Malignant_astrolike","Malignant_cycling","Malignant_intermediate",      
              "Malignant_OPClike","Malignant_reactive","Myeloid_Inflammatory","Myeloid_Macrophage",
              "Myeloid_MacrophageSPP1","Myeloid_Microglia","Neuron")
results <- lapply(features, function(x) {
  fit <- lm(reformulate(c("Grade", "treatment", "ini_rec"), response = x), data = lm_df.oligo)
  
  co <- summary(fit)$coefficients
  ci <- confint(fit)
  
  data.frame(feature = x, grade_beta = co["Grade", "Estimate"], grade_p = co["Grade", "Pr(>|t|)"], grade_low = ci["Grade", 1],
             grade_high = ci["Grade", 2], treatmentPCV_beta = co["treatmentPCV", "Estimate"],
             treatmentPCV_p = co["treatmentPCV", "Pr(>|t|)"], treatmentPCV_low = ci["treatmentPCV", 1],
             treatmentPCV_high = ci["treatmentPCV", 2], treatmentTMZ_beta = co["treatmentTMZ", "Estimate"],
             treatmentTMZ_p = co["treatmentTMZ", "Pr(>|t|)"], treatmentTMZ_low = ci["treatmentTMZ", 1],
             treatmentTMZ_high = ci["treatmentTMZ", 2])
})
results <- do.call(rbind, results)
plot_df <- results %>%
  select(feature, grade_beta, grade_p, grade_low, grade_high,
    treatmentPCV_beta, treatmentPCV_p, treatmentPCV_low, treatmentPCV_high,
    treatmentTMZ_beta, treatmentTMZ_p, treatmentTMZ_low, treatmentTMZ_high) %>%
  pivot_longer(cols = -feature, names_to = c("variable", ".value"),
    names_pattern = "(grade|treatmentPCV|treatmentTMZ)_(.*)")
plot_df$feature = gsub("like","-like",plot_df$feature)
plot_df$feature = gsub("SPP1","SPP1+",plot_df$feature)
plot_df$feature = factor(plot_df$feature, 
                              levels = c(unique(plot_df$feature)[!grepl("^M",unique(plot_df$feature))], 
                                         unique(plot_df$feature)[grepl("Myeloid_",unique(plot_df$feature))],
                                         unique(plot_df$feature)[grepl("Malignant_",unique(plot_df$feature))]))
#plot figure of coefficient estimates and confidence intervals
ggplot(plot_df, aes(x = beta, y = feature, color = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar( aes(xmin = low, xmax = high), orientation = "y", height = 0,
                 position = position_dodge(width = 0.35)) +
  geom_point( size = 3, position = position_dodge(width = 0.35)) +
  theme_classic() + labs( x = "Regression coefficient", y = NULL, color = NULL)
ggsave(paste0(fig.dir,"/lm_estimates_oligo_full_clusters.pdf"), width = 6, height = 5)



##barplot ordered by PC1 and 2

#make colour scales
cols.mut = c("1" = "grey30", "0" = "gray95")
cols.grade = c("2" = "skyblue", "3" = "royalblue", "4" = "darkblue")
cols.ini_rec = c("ini" = "violet", "rec1" = "purple", "rec2" = "purple4")
full.color.scale = c(cols.mut, cols.grade, cols.ini_rec)

#make annotation dataframe
annot_data = data.frame(sample_id = unique(df_counts$sample))
annot_data$sample = sapply(annot_data$sample_id, function(x) {
  samp.md$Sample[which(samp.md$sample_id==x)]
})
annot_data$cancer = sapply(annot_data$sample, function(x) {unlist(strsplit(x,"_"))[1]})
annot_data$ini_rec = sapply(annot_data$sample, function(x) {unlist(strsplit(x,"_"))[3]})
annot_data$grade = sapply(annot_data$sample, function(x) {
  samp.md$Grade[which(samp.md$Sample==x)]
})
annot_data$PCV = sapply(annot_data$sample, function(x) {
  samp.md$PCV[which(samp.md$Sample==x)]
})
annot_data$TMZ = sapply(annot_data$sample, function(x) {
  samp.md$TMZ[which(samp.md$Sample==x)]
})
annot_data$Radiotherapy = sapply(annot_data$sample, function(x) {
  samp.md$Radiotherapy[which(samp.md$Sample==x)]
})
annot_data$sample = annot_data$sample_id


#make survival time dataframe
rec_time_data = data.frame(sample_id = unique(df_counts$sample))
rec_time_data$sample = sapply(rec_time_data$sample_id, function(x) {
  samp.md$Sample[which(samp.md$sample_id==x)]
})
rec_time_data$cancer = sapply(rec_time_data$sample, function(x) {unlist(strsplit(x,"_"))[1]})
rec_time_data$patient = sapply(rec_time_data$sample, function(x) {unlist(strsplit(x,"_"))[2]})
rec_time_data$ini_rec = sapply(rec_time_data$sample, function(x) {unlist(strsplit(x,"_"))[3]})

rec_time_data$surv_time = sapply(seq(length(rec_time_data$patient)), function(x) {
  rel.rows = samp.md$DDN[which(samp.md$patient==rec_time_data$patient[x])]
  as.numeric(rel.rows[which(!is.na(rel.rows))] - samp.md$Surgery_date[which(samp.md$patient==rec_time_data$patient[x] & samp.md$ini_rec==rec_time_data$ini_rec[x])])/365
})

#make barplot of astro samples ordered by PC1 value
df_to_plot = df_counts[df_counts$cancer=="Astrocytoma",]
df_to_plot <- df_to_plot %>% group_by(sample) %>%
  mutate(count = ifelse(cell.type == "Myeloid_Macrophage",
      count + sum(count[cell.type == "Myeloid_Macrophage SPP1+"]), count)) %>%
  filter(cell.type != "Myeloid_Macrophage SPP1+")
samples_order = unique(df_to_plot$sample)[order(-pca.astro$x[,1])]
df_to_plot$sample = factor(df_to_plot$sample, levels = samples_order)
cell.type.pal = gg_color_hue(length(unique(df_to_plot$cell.type)))
names(cell.type.pal) = levels(df_to_plot$cell.type)
annot_to_plot = annot_data[annot_data$cancer=="astro",]
annot_to_plot$sample = factor(annot_to_plot$sample, levels = samples_order)
feats.keep = c("ini_rec", "grade", "TMZ", "PCV", "Radiotherapy")
annot_long <- annot_to_plot %>%
  mutate(across(-sample, as.character)) %>%   # convert all annotation vars to character
  pivot_longer(cols = -sample, names_to = "feature", values_to = "value") %>%
  filter(feature %in% feats.keep) %>%   # keep only desired features
  mutate(sample = factor(sample, levels = samples_order),
         feature = factor(feature, levels = rev(feats.keep)),
         value  = factor(value) )

rec_to_plot = rec_time_data[which(rec_time_data$cancer=="astro"),]
rec_to_plot$sample = rec_to_plot$sample_id
rec_to_plot$sample = factor(rec_to_plot$sample, levels = samples_order)
df_to_plot$cell.type[which(df_to_plot$cell.type=="Malignant_reactive")] = "Malignant_reactive/MES-like"
df_to_plot$cell.type = factor(df_to_plot$cell.type, levels = rev(c("Malignant_reactive/MES-like","Malignant_astro-like","Malignant_intermediate","Malignant_neuro-like","Malignant_OPC-like","Malignant_cycling",
                                                               "Myeloid_Inflammatory","Myeloid_Microglia","Myeloid_Macrophage",
                                                               "Neuron","Glial")))
state.cols = gg_color_hue(4)
names(state.cols) = c("Malignant_astro-like","Malignant_intermediate","Malignant_neuro-like","Malignant_OPC-like")
myeloid.cols = c("#8B1A1A","#F4A582","#4D7BA5")
names(myeloid.cols) = c("Myeloid_Inflammatory","Myeloid_Microglia","Myeloid_Macrophage")
other.malig.cols = c("grey50","gold")
names(other.malig.cols) = c("Malignant_cycling","Malignant_reactive/MES-like")
tme.cols = c("darkseagreen","plum1")
names(tme.cols) = c("Glial","Neuron")
cell.type.cols = c(state.cols,myeloid.cols,other.malig.cols,tme.cols)


p1 = ggplot(annot_long, aes(x = sample, y = feature, fill = value)) +
  geom_tile(color = "black") + scale_fill_manual(values = full.color.scale) +
  labs(x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),legend.position = "none", text=element_text(size=14))

p2 = ggplot(rec_to_plot, aes(x = sample, y = surv_time)) +
  geom_linerange(aes(ymin = 0, ymax = surv_time)) + geom_point(size=2) + labs(x="",y = "Survival time")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text=element_text(size=14))

p3 = ggplot(df_to_plot, aes(fill=cell.type, y=count, x=sample)) + scale_fill_manual(values=cell.type.pal)+
  geom_bar(position="fill", stat="identity")+ labs(fill='Cell type',x = "sample", y = "Proportion") +
  theme(text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = cell.type.cols)

wrap_plots(list(p1,p2,p3), ncol=1, heights = c(1,2,3))
ggsave(paste0(fig.dir,"/barplot_celltypes_by_PC1_astro_full_clusters.pdf"), height = 5, width = 10)


samples_order = unique(df_to_plot$sample)[order(pca.astro$x[,2])]
df_to_plot$sample = factor(df_to_plot$sample, levels = samples_order)

annot_to_plot = annot_data[annot_data$cancer=="astro",]
annot_to_plot$sample = factor(annot_to_plot$sample, levels = samples_order)
feats.keep = c("ini_rec", "grade", "TMZ", "PCV", "Radiotherapy")
annot_long <- annot_to_plot %>%
  mutate(across(-sample, as.character)) %>%   # convert all annotation vars to character
  pivot_longer(cols = -sample,names_to = "feature",values_to = "value") %>%
  filter(feature %in% feats.keep) %>%   # keep only desired features
  mutate(sample = factor(sample, levels = samples_order),
         feature = factor(feature, levels = rev(feats.keep)),
         value  = factor(value))

rec_to_plot = rec_time_data[which(rec_time_data$cancer=="astro"),]
rec_to_plot$sample = rec_to_plot$sample_id
rec_to_plot$sample = factor(rec_to_plot$sample, levels = samples_order)
df_to_plot$cell.type[which(df_to_plot$cell.type=="Malignant_reactive")] = "Malignant_reactive/MES-like"
df_to_plot$cell.type = factor(df_to_plot$cell.type, levels = rev(c("Malignant_reactive/MES-like","Malignant_astro-like","Malignant_intermediate","Malignant_neuro-like","Malignant_OPC-like","Malignant_cycling",
                                                                   "Myeloid_Inflammatory","Myeloid_Microglia","Myeloid_Macrophage",
                                                                   "Neuron","Glial")))
p1 = ggplot(annot_long, aes(x = sample, y = feature, fill = value)) +
  geom_tile(color = "black") + scale_fill_manual(values = full.color.scale) + labs(x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),legend.position = "none", text=element_text(size=14))

p2 = ggplot(rec_to_plot, aes(x = sample, y = surv_time)) +
  geom_linerange(aes(ymin = 0, ymax = surv_time)) + geom_point(size=2) + labs(x="",y = "Survival time")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text=element_text(size=14))

p3 = ggplot(df_to_plot, aes(fill=cell.type, y=count, x=sample)) +scale_fill_manual(values=cell.type.pal) +
  geom_bar(position="fill", stat="identity")+ labs(fill='Cell type',x = "sample", y = "Proportion") +
  theme(text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = cell.type.cols)

wrap_plots(list(p1,p2,p3), ncol=1, heights = c(1,1,3))
ggsave(paste0(fig.dir,"/barplot_celltypes_by_PC2_astro_full_clusters.pdf"), height = 5, width = 10)





#same for oligos
df_to_plot = df_counts[df_counts$cancer=="Oligodendroglioma",]
df_to_plot <- df_to_plot %>% group_by(sample) %>%
  mutate(count = ifelse(cell.type == "Myeloid_Macrophage",
                        count + sum(count[cell.type == "Myeloid_Macrophage SPP1+"]), count)) %>%
  filter(cell.type != "Myeloid_Macrophage SPP1+")
samples_order = unique(df_to_plot$sample)[order(pca.oligo$x[,1])]
df_to_plot$sample = factor(df_to_plot$sample, levels = samples_order)

annot_to_plot = annot_data[annot_data$cancer=="oligo",]
annot_to_plot$sample = factor(annot_to_plot$sample, levels = samples_order)
feats.keep = c("ini_rec", "grade", "TMZ", "PCV", "Radiotherapy")
annot_long <- annot_to_plot %>%
  mutate(across(-sample, as.character)) %>%   # convert all annotation vars to character
  pivot_longer(cols = -sample,names_to = "feature",values_to = "value") %>%
  filter(feature %in% feats.keep) %>%   # keep only desired features
  mutate(sample = factor(sample, levels = samples_order),
         feature = factor(feature, levels = rev(feats.keep)),value  = factor(value))

rec_to_plot = rec_time_data[which(rec_time_data$cancer=="oligo"),]
rec_to_plot$state = sapply(seq(length(rec_to_plot$patient)), function(x) {
  rel.rows = samp.md$State_DDN[which(samp.md$patient==rec_to_plot$patient[x])]
  state = rel.rows[!is.na(rel.rows)]
  if(state=="Dead"){
    "Dead"
  }else{
    "Alive"
  }
})
rec_to_plot$state = factor(rec_to_plot$state, levels = c("Alive","Dead"))
rec_to_plot$sample = rec_to_plot$sample_id
rec_to_plot$sample = factor(rec_to_plot$sample, levels = samples_order)
df_to_plot$cell.type[which(df_to_plot$cell.type=="Malignant_reactive")] = "Malignant_reactive/MES-like"
df_to_plot$cell.type = factor(df_to_plot$cell.type, levels = rev(c("Malignant_reactive/MES-like","Malignant_astro-like","Malignant_intermediate","Malignant_OPC-like","Malignant_cycling",
                                                                   "Myeloid_Inflammatory","Myeloid_Microglia","Myeloid_Macrophage",
                                                                   "Neuron","Glial")))
p1 = ggplot(annot_long, aes(x = sample, y = feature, fill = value)) +
  geom_tile(color = "black") + scale_fill_manual(values = full.color.scale) +
  labs(x = "", y = "") + theme_minimal() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),legend.position = "none", text=element_text(size=14))

p2 = ggplot(rec_to_plot, aes(x = sample, y = surv_time)) +
  geom_linerange(aes(ymin = 0, ymax = surv_time)) + geom_point(size=1.5, aes(color=state)) + labs(x="",y = "Survival time")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text=element_text(size=14))

p3 = ggplot(df_to_plot, aes(fill=cell.type, y=count, x=sample)) +scale_fill_manual(values=cell.type.pal)+
  geom_bar(position="fill", stat="identity")+ labs(fill='Cell type',x = "sample", y = "Proportion") +
  theme(text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = cell.type.cols)

wrap_plots(list(p1,p2,p3), ncol=1, heights = c(1,2,3))
ggsave(paste0(fig.dir,"/barplot_celltypes_by_PC1_oligo_full_clusters.pdf"), height = 5, width = 10)


