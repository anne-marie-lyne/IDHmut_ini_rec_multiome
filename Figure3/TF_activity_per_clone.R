library(Seurat)
library(Signac)
library(infercnv, lib.loc="/mnt/beegfs/home/alyne/persistent/R_libs")
library(data.table)
library(GenomicRanges)
library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)
library(miscTools)
library(ggpubr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


out.dir = "/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_oligo/figures_scenicplus_newdb_oligo"
fig.dir = "/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_oligo/figures_TF_act_by_clone"
dir.create(fig.dir)

#first get chromosome of all TF in scenicplus output
#read in gene AUC scenicplus output
use.extended = F
TF.AUCg.d = fread(paste0(out.dir,"/direct_gene_based_AUC_oligo.tsv.gz"), header=T, data.table = F)
if(use.extended){
  TF.AUCg.e = fread(paste0(out.dir,"/extended_gene_based_AUC_oligo.tsv.gz"), header=T, data.table = F)
}else{
  TF.AUCg.e = data.frame(Cell = TF.AUCg.d$Cell)
}
TF.AUCg = cbind(TF.AUCg.d,TF.AUCg.e[,-1])
rownames(TF.AUCg) = TF.AUCg$Cell
rownames(TF.AUCg) = sapply(rownames(TF.AUCg), function(x) {unlist(strsplit(x,"___",fixed=T))[1]})
TF.AUCg$Cell = NULL

#for now keep only +/+
TF.AUCg = TF.AUCg[,grepl("\\+/\\+",colnames(TF.AUCg))]

TFs.all = unique(sapply(colnames(TF.AUCg), function(x) {unlist(strsplit(x,"_"))[1]}))


#read in gene info
gene.info = fread("/mnt/beegfs/common/apps/cellranger/refdata/refdata-gex-GRCh38_and_GRCm39-2024-A/genes/genes.gtf.gz", 
                  data.table=F, header=F)
gene.info = gene.info[gene.info$V3=="gene",]
gene.info = gene.info[grepl("protein_coding",gene.info$V9),]

rel.rows = lapply(TFs.all, function(x) {gene.info[grepl(paste0("\"GRCh38_",x,"\""),gene.info$V9),]})
TFs.chr = sapply(rel.rows, function(x) {
  if(nrow(x)>0){
    unlist(strsplit(x[,1],"_"))[2]
  }else{
    NA
  }
})

# TFs.poss = TFs.all[which(TFs.chr %in% c("chr4","chr7","chr17"))]
TFs.poss = TFs.all
TF.inds = sort(unique(unlist(lapply(TFs.poss, 
                                    function(x) {which(grepl(paste0("^",x), colnames(TF.AUCg)))}))))

#for each clone/sample, plot activity/accessibility and expression by cell state
CNVs <- as.data.frame(read_excel("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/summary_main_CNA_all_samples.xlsx", na = "NA"))
CNVs.oligo = CNVs[grepl("oligo",CNVs$sample),]



#plot figure to compare to infercnv
CNV.props.mat <- as.data.frame(read_excel("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/cell_states_and_CNA_per_clone.xlsx", 
                                          col_types = "text", na = "NA"))
CNV.props.mat[is.na(CNV.props.mat)] = 0
CNV.props.mat[-1] <- lapply(CNV.props.mat[-1], as.numeric)

cell.states = c("astro-like","intermediate","OPC-like")
col.states = gg_color_hue(length(cell.states))
names(col.states) = cell.states

gene.ord = fread(paste0("/mnt/beegfs/home/alyne/stagein/data_tmp/numbat_inputs/infercnv_atac_oligo",
                        "/gene_order_oligo.tsv"), header=F, data.table = F, fill=TRUE)



#in first approach, ignore cell states
df.clone = list()
oligo.list = list()
icnv.list = list()
count = 1
for(i in 1:length(CNVs.oligo$sample)){
# for(i in 1:6){
  print(i)
  samp = CNVs.oligo$sample[i]
  oligo.list[[i]] = readRDS(paste0("/mnt/beegfs/home/alyne/stagein/oligo_seurat/oligo_malig_",samp,".rds"))
  
  if(!is.na(CNVs.oligo$CNA_use_clustering[which(CNVs.oligo$sample==samp)])){
    #sample has clones, so do clustering
    out.dir = paste0("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/infercnv_ATAC_reference_",samp)
    if(file.exists(paste0(out.dir,"/run.final.infercnv_obj"))){
      obj <- readRDS(paste0(out.dir,"/run.final.infercnv_obj"))
      obs = obj@expr.data
      to.plot = as.matrix(t(obs)[unlist(obj@observation_grouped_cell_indices),])
      
      #find which chromosomes to use
      chr.to.plot = paste0("chr",unlist(strsplit(CNVs.oligo$CNA_use_clustering[which(CNVs.oligo$sample==samp)],";")))
      
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
      icnv.list[[count]] <- CreateSeuratObject(assay.v5)
      # icnv <- FindVariableFeatures(icnv, selection.method = "vst", nfeatures = 3000)
      all.genes <- rownames(icnv.list[[count]])
      icnv.list[[count]] <- ScaleData(icnv.list[[count]], features = all.genes)
      icnv.list[[count]] = RunPCA(icnv.list[[count]], features = all.genes)
      icnv.list[[count]] <- FindNeighbors(icnv.list[[count]], dims = 1:10)
      icnv.list[[count]] <- FindClusters(icnv.list[[count]], resolution = 0.1)
      if(samp %in% c("oligo_P0007_ini_S0015","astro_P0026_ini_S0050","oligo_P0031_ini_S0060")){
        icnv.list[[count]] <- FindClusters(icnv.list[[count]], resolution = 0.2)
      }
      if(samp %in% c("oligo_P0005_rec1_S0012")){
        icnv.list[[count]] <- FindClusters(icnv.list[[count]], resolution = 0.3)
      }
      icnv.list[[count]] <- RunUMAP(icnv.list[[count]], dims = 1:10)
      ind = match(colnames(icnv.list[[count]]),colnames(oligo.list[[i]]))
      icnv.list[[count]]$ICA_score_states = oligo.list[[i]]$ICA_score_states[ind]
      
      df.clone.sub = list()
      #for each clone, make data frame summarising states and TF activity/expression
      for(j in 1:length(levels(icnv.list[[count]]$seurat_clusters))){
        print(paste("clone",j))
        df.samp.info = data.frame(sample = paste0(samp,"_",levels(icnv.list[[count]]$seurat_clusters)[j]),
                                  Nastro = sum(icnv.list[[count]]$ICA_score_states=="astro-like" & icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j]),
                                  Nopc = sum(icnv.list[[count]]$ICA_score_states=="OPC-like" & icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j]),
                                  Nint = sum(icnv.list[[count]]$ICA_score_states=="intermediate" & icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j]),
                                  Ntot = sum(icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j]))
        #get summary of TF expression
        row.inds = match(TFs.poss, rownames(oligo.list[[i]]@assays$SCT@data))
        cells_clone = colnames(icnv.list[[count]])[which(icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j])]
        col.inds = match(cells_clone, colnames(oligo.list[[i]]))
        expr_means = sapply(row.inds, function(x) {
          if(!is.na(x)){
            mean(oligo.list[[i]]@assays$SCT@data[x,col.inds])
          }else{NA}
        })
        names(expr_means) = paste0(TFs.poss,"_expr_mean")
        expr_meds = sapply(row.inds, function(x) {
          if(!is.na(x)){
            median(oligo.list[[i]]@assays$SCT@data[x,col.inds])
          }else{NA}
        })
        names(expr_meds) = paste0(TFs.poss,"_expr_median")
        #get summary of TF activity
        row.inds = match(paste0(cells_clone,"-",samp),rownames(TF.AUCg))
        act_means = colMeans(TF.AUCg[row.inds,TF.inds], na.rm=T)
        names(act_means) = paste0(names(act_means),"_act_mean")
        act_meds = colMedians(TF.AUCg[row.inds,TF.inds], na.rm=T)
        names(act_meds) = paste0(names(act_meds),"_act_median")
        df.TF = data.frame(t(c(expr_means,expr_meds,act_means,act_meds)))
        #get summary of CNA
        ind.cnv = which(CNV.props.mat$Clone==paste0(samp,"_",levels(icnv.list[[count]]$seurat_clusters)[j]))
        CNV.name = paste(colnames(CNV.props.mat[,6:ncol(CNV.props.mat)])[which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
        CNV.direc = paste(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)][which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
        df.CNV = data.frame(CNV.name = CNV.name, CNV.direc = CNV.direc, 
                            chr4 = grepl("chr4",CNV.name), chr7 = grepl("chr7",CNV.name), chr11 = grepl("chr11",CNV.name),
                            chr13 = grepl("chr13",CNV.name), chr14 = grepl("chr14",CNV.name),
                            chr17 = grepl("chr17",CNV.name), chr18 = grepl("chr18",CNV.name))
        df.clone.sub[[j]] = cbind(df.samp.info,df.TF,df.CNV)
      }
      df.clone[[i]] = do.call("rbind",df.clone.sub)
      count = count + 1
    }else{
      df.clone[[i]] = NA
    }
  }else{
    #no clones
    #get summary of cell state counts
    print("got to here")
    df.samp.info = data.frame(sample = samp,Nastro = sum(oligo.list[[i]]$ICA_score_states=="astro-like"),
                              Nopc = sum(oligo.list[[i]]$ICA_score_states=="OPC-like"),
                              Nint = sum(oligo.list[[i]]$ICA_score_states=="intermediate"), Ntot = ncol(oligo.list[[i]]))
    #get summary of TF expression
    row.inds = match(TFs.poss, rownames(oligo.list[[i]]@assays$SCT@data))
    print("now to here")
    expr_means = sapply(row.inds, function(x) {
      if(!is.na(x)){
        mean(oligo.list[[i]]@assays$SCT@data[x,])
      }else{NA}
    })
    names(expr_means) = paste0(TFs.poss,"_expr_mean")
    expr_meds = sapply(row.inds, function(x) {
      if(!is.na(x)){
        median(oligo.list[[i]]@assays$SCT@data[x,])
      }else{NA}
    })
    print("next to here")
    names(expr_meds) = paste0(TFs.poss,"_expr_median")
    #get summary of TF activity
    row.inds = match(paste0(colnames(oligo.list[[i]]),"-",samp),rownames(TF.AUCg))
    act_means = colMeans(TF.AUCg[row.inds,TF.inds],na.rm = T)
    names(act_means) = paste0(names(act_means),"_act_mean")
    act_meds = colMedians(TF.AUCg[row.inds,TF.inds],na.rm = T)
    names(act_meds) = paste0(names(act_meds),"_act_median")
    df.TF = data.frame(t(c(expr_means,expr_meds,act_means,act_meds)))
    print("nearly finished")
    #get summary of CNA
    ind.cnv = which(CNV.props.mat$Clone==samp)
    CNV.name = paste(colnames(CNV.props.mat[,6:ncol(CNV.props.mat)])[which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
    CNV.direc = paste(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)][which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
    df.CNV = data.frame(CNV.name = CNV.name, CNV.direc = CNV.direc, 
                        chr4 = grepl("chr4",CNV.name), chr7 = grepl("chr7",CNV.name), chr11 = grepl("chr11",CNV.name),
                        chr13 = grepl("chr13",CNV.name), chr14 = grepl("chr14",CNV.name),
                        chr17 = grepl("chr17",CNV.name), chr18 = grepl("chr18",CNV.name))
    df.clone[[i]] = cbind(df.samp.info,df.TF,df.CNV)
    print("finished")
  }
}
df.clone.all = do.call("rbind",df.clone)

fwrite(df.clone.all, paste0(fig.dir,"/TF_activity_by_clone_oligos.tsv.gz"), 
       sep ="\t", row.names = T, col.names = T)

df.clone.all = fread(paste0(fig.dir,"/TF_activity_by_clone_oligos.tsv.gz"), header=T, data.table=F)


#same thing but splitting cells into malignant states
#in first approach, ignore cell states
df.clone = list()
count = 1
for(i in 1:length(CNVs.oligo$sample)){
  # for(i in 1:6){
  print(i)
  samp = CNVs.oligo$sample[i]
  
  if(!is.na(CNVs.oligo$CNA_use_clustering[which(CNVs.oligo$sample==samp)])){
    #sample has clones, so do clustering
    out.dir = paste0("/mnt/beegfs/home/alyne/stageout/infercnv_outputs/infercnv_ATAC_reference_",samp)
    if(file.exists(paste0(out.dir,"/run.final.infercnv_obj"))){
      df.clone.sub = list()
      #for each clone, make data frame summarising states and TF activity/expression
      for(j in 1:length(levels(icnv.list[[count]]$seurat_clusters))){
        print(paste("clone",j))
        df.samp.info = data.frame(sample = paste0(samp,"_",levels(icnv.list[[count]]$seurat_clusters)[j]))
        #get summary of TF expression
        df.clone.sub.state = list()
        for(k in 1:3){
          cells_clone = colnames(icnv.list[[count]])[which(icnv.list[[count]]$seurat_clusters==levels(icnv.list[[count]]$seurat_clusters)[j]) & 
                                                       icnv.list[[count]]$ICA_score_states==cell.states[k]]
          row.inds = match(TFs.poss, rownames(oligo.list[[i]]@assays$SCT@data))
          col.inds = match(cells_clone, colnames(oligo.list[[i]]))
          expr_means = sapply(row.inds, function(x) {
            if(!is.na(x)){
              mean(oligo.list[[i]]@assays$SCT@data[x,col.inds])
            }else{NA}
          })
          names(expr_means) = paste0(TFs.poss,"_expr_mean")
          expr_meds = sapply(row.inds, function(x) {
            if(!is.na(x)){
              median(oligo.list[[i]]@assays$SCT@data[x,col.inds])
            }else{NA}
          })
          names(expr_meds) = paste0(TFs.poss,"_expr_median")
          #get summary of TF activity
          row.inds = match(paste0(cells_clone,"-",samp),rownames(TF.AUCg))
          act_means = colMeans(TF.AUCg[row.inds,TF.inds], na.rm=T)
          names(act_means) = paste0(names(act_means),"_act_mean")
          act_meds = colMedians(TF.AUCg[row.inds,TF.inds], na.rm=T)
          names(act_meds) = paste0(names(act_meds),"_act_median")
          df.TF = data.frame(t(c(expr_means,expr_meds,act_means,act_meds)))
          #get summary of CNA
          ind.cnv = which(CNV.props.mat$Clone==paste0(samp,"_",levels(icnv.list[[count]]$seurat_clusters)[j]))
          CNV.name = paste(colnames(CNV.props.mat[,6:ncol(CNV.props.mat)])[which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
          CNV.direc = paste(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)][which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
          df.CNV = data.frame(CNV.name = CNV.name, CNV.direc = CNV.direc, state = cell.states[k],
                              chr4 = grepl("chr4",CNV.name), chr7 = grepl("chr7",CNV.name), chr11 = grepl("chr11",CNV.name),
                              chr13 = grepl("chr13",CNV.name), chr14 = grepl("chr14",CNV.name),
                              chr17 = grepl("chr17",CNV.name), chr18 = grepl("chr18",CNV.name))
          df.clone.sub.state[[k]] = cbind(df.samp.info,df.TF,df.CNV)
        }
        df.clone.sub[[j]] = do.call("rbind",df.clone.sub.state)
      }
      df.clone[[i]] = do.call("rbind",df.clone.sub)
      count = count + 1
    }else{
      df.clone[[i]] = NA
    }
  }else{
    #no clones
    df.samp.info = data.frame(sample = samp)
    #get summary of TF expression
    df.clone.state = list()
    for(k in 1:3){
      #get cells in this state
      row.inds = match(TFs.poss, rownames(oligo.list[[i]]@assays$SCT@data))
      col.inds = which(oligo.list[[i]]$ICA_score_states==cell.states[k])
      expr_means = sapply(row.inds, function(x) {
        if(!is.na(x)){
          mean(oligo.list[[i]]@assays$SCT@data[x,col.inds])
        }else{NA}
      })
      names(expr_means) = paste0(TFs.poss,"_expr_mean")
      expr_meds = sapply(row.inds, function(x) {
        if(!is.na(x)){
          median(oligo.list[[i]]@assays$SCT@data[x,col.inds])
        }else{NA}
      })
      names(expr_meds) = paste0(TFs.poss,"_expr_median")
      #get summary of TF activity
      row.inds = match(paste0(colnames(oligo.list[[i]])[col.inds],"-",samp),rownames(TF.AUCg))
      act_means = colMeans(TF.AUCg[row.inds,TF.inds],na.rm = T)
      names(act_means) = paste0(names(act_means),"_act_mean")
      act_meds = colMedians(TF.AUCg[row.inds,TF.inds],na.rm = T)
      names(act_meds) = paste0(names(act_meds),"_act_median")
      df.TF = data.frame(t(c(expr_means,expr_meds,act_means,act_meds)))
      #get summary of CNA
      ind.cnv = which(CNV.props.mat$Clone==samp)
      CNV.name = paste(colnames(CNV.props.mat[,6:ncol(CNV.props.mat)])[which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
      CNV.direc = paste(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)][which(abs(CNV.props.mat[ind.cnv,6:ncol(CNV.props.mat)])==1)],collapse=";")
      df.CNV = data.frame(CNV.name = CNV.name, CNV.direc = CNV.direc, state = cell.states[k],
                          chr4 = grepl("chr4",CNV.name), chr7 = grepl("chr7",CNV.name), chr11 = grepl("chr11",CNV.name),
                          chr13 = grepl("chr13",CNV.name), chr14 = grepl("chr14",CNV.name),
                          chr17 = grepl("chr17",CNV.name), chr18 = grepl("chr18",CNV.name))
      df.clone.state[[k]] = cbind(df.samp.info,df.TF,df.CNV)
    }
    df.clone[[i]] = do.call("rbind",df.clone.state)
  }
}
df.clone.all.bystate = do.call("rbind",df.clone)
fwrite(df.clone.all.bystate, paste0(fig.dir,"/TF_activity_by_state_by_clone_oligos.tsv.gz"), 
       sep ="\t", row.names = T, col.names = T)


df.clone.all.bystate = fread(paste0(fig.dir,"/TF_activity_by_state_by_clone_oligos.tsv.gz"), header=T, data.table=F)
df.clone.all = fread(paste0(fig.dir,"/TF_activity_by_clone_oligos.tsv.gz"), header=T, data.table=F)

#add Ntot
df.clone.all.bystate$Ntot = unlist(lapply(df.clone.all$Ntot, function(x) {rep(x,3)}))




#compute p values using linear model
library(lme4)
library(lmerTest)

# linear model to include sample effect
# but won't solve for TF expression so only include sample effect for activity
Ncell.min=100
lmer.expr = list()
for(x in seq(length(TFs.all))){
  act.colnames = colnames(df.clone.all.bystate)[grepl(paste0("^",TFs.all[x],"_direct"),colnames(df.clone.all.bystate))]
  if(length(act.colnames)>0){
    df.to.plot = df.clone.all.bystate[df.clone.all.bystate$Ntot>=Ncell.min,c("sample","chr4",
                                                                             paste0(TFs.all[x],"_expr_mean"),act.colnames[1],"state")]
    df.to.plot = df.to.plot[df.to.plot$state=="intermediate",]
    colnames(df.to.plot) = c("clone","chromosome","expression_mean","activity_mean","state")
    df.to.plot$sample = sapply(df.to.plot$clone, function(x) {unlist(strsplit(x,"_"))[4]})
    #keep only full samples
    
    df.to.plot = df.to.plot[!is.na(df.to.plot$expression_mean), ]
    lmer.expr[[x]] = tryCatch(
      lm(expression_mean ~ chromosome, data = df.to.plot),
      error = function(e) NA
    )
  }else{lmer.expr[[x]] = NA}
}
names(lmer.expr) = TFs.all

lmer.act = list()
for(x in seq(length(TFs.all))){
  act.colnames = colnames(df.clone.all)[grepl(paste0("^",TFs.all[x],"_direct"),colnames(df.clone.all))]
  if(length(act.colnames)>0){
    df.to.plot = df.clone.all[df.clone.all$Ntot>=Ncell.min,c("sample","chr4",paste0(TFs.all[x],"_expr_mean"),act.colnames[1])]
    colnames(df.to.plot) = c("clone","chromosome","expression_mean","activity_mean")
    df.to.plot$sample = sapply(df.to.plot$clone, function(x) {unlist(strsplit(x,"_"))[4]})
    df.to.plot = df.to.plot[!is.na(df.to.plot$activity_mean), ]
    lmer.act[[x]] = lmerTest::lmer(activity_mean ~ chromosome + (1 | sample), data = df.to.plot)
  }else{lmer.act[[x]] = NA}
}
names(lmer.act) = TFs.all


#extract p-values and BH correct
pval.expr <- sapply(lmer.expr, function(x) {
  if (inherits(x, "lm")) {
    summary(x)$coefficients[2,4]
  }else{ NA_real_ }
})
pval.expr.adj <- p.adjust(pval.expr, method = "BH")
pval.act <- sapply(lmer.act, function(x) {
  if (inherits(x, "lmerMod")) {
    anova(x)["chromosome", "Pr(>F)"]
  }else{ NA_real_ }
})
pval.act.adj <- p.adjust(pval.act, method = "BH")

#compute "mean" p-value across bth analyses for each TF
#try to identify those with diffrent expression and activity
mean.pvals = sapply(1:length(pval.act.adj), function(x) {if(any(is.na(c(pval.expr.adj[x],pval.act.adj[x])))){
  NA
}else{
  mean(-1*log10(c(pval.expr.adj[x],pval.act.adj[x])))
}
})
names(mean.pvals) = TFs.all
mean.pvals = sort(mean.pvals,decreasing=T,na.last=T)



#make final plot, remove anything with expression and activity in different directions
#include those from ~top 25 on chr4, astrocyte-lineage, OPC lineage, and cycling
TFs.to.plot = c("SMAD1","PDLIM5","NR3C2","SOX9","SOX6","ETV1","ETS1","E2F3")

df.list = list()
for(i in 1:length(TFs.to.plot)){
  act.colnames = colnames(df.clone.all.bystate)[grepl(paste0("^",TFs.to.plot[i],"_direct"),colnames(df.clone.all.bystate))]
  df.list[[i]] = df.clone.all.bystate[df.clone.all.bystate$Ntot>=Ncell.min,c("sample","chr4",act.colnames[1],"state")]
  df.list[[i]] = df.list[[i]][df.list[[i]]$state=="intermediate",]
  df.list[[i]]$TF = TFs.to.plot[i]
  colnames(df.list[[i]]) = c("sample","chr4","activity_mean","state","TF")
}
df.act = do.call("rbind",df.list)
df.act$TF = factor(df.act$TF, levels = TFs.to.plot)
df.act$chr4[which(df.act$chr4=="TRUE")] = "deletion"
df.act$chr4[which(df.act$chr4=="FALSE")] = "no deletion"
df.act$chr4 = factor(df.act$chr4, levels = c("no deletion","deletion"))

df.list2 = list()
for(i in 1:length(TFs.to.plot)){
  df.list2[[i]] = df.clone.all.bystate[df.clone.all.bystate$Ntot>=Ncell.min,c("sample","chr4",paste0(TFs.to.plot[i],"_expr_mean"),"state")]
  df.list2[[i]] = df.list2[[i]][df.list2[[i]]$state=="intermediate",]
  df.list2[[i]]$TF = TFs.to.plot[i]
  colnames(df.list2[[i]]) = c("sample","chr4","expression_mean","state","TF")
}
df.expr = do.call("rbind",df.list2)
df.expr$TF = factor(df.expr$TF, levels = TFs.to.plot)
df.expr$chr4[which(df.expr$chr4=="TRUE")] = "deletion"
df.expr$chr4[which(df.expr$chr4=="FALSE")] = "no deletion"
df.expr$chr4 = factor(df.expr$chr4, levels = c("no deletion","deletion"))

pvals.expr = data.frame(TF = TFs.to.plot, pvalue = paste(signif(pval.expr.adj[TFs.to.plot],2)), y = max(df.expr$expression_mean)+0.2)
pvals.act = data.frame(TF = TFs.to.plot, pvalue = paste(signif(pval.act.adj[TFs.to.plot],2)), y = max(df.act$activity_mean,na.rm=T)+0.01)

p1 = ggplot(df.expr, aes(x = TF, y = expression_mean, fill = chr4)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.08,dodge.width = 0.75),size = 1) +
  labs(y = "Mean expression",x="") +
  geom_text(data = pvals.expr,aes(x = TF, y = y, label = pvalue),inherit.aes = FALSE,size = 3)
p2 = ggplot(df.act, aes(x = TF, y = activity_mean, fill = chr4)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.08,dodge.width = 0.75),size = 1) +
  labs(y = "Mean activity",x="") + guides(fill="none") +
  geom_text(data = pvals.act,aes(x = TF, y = y, label = pvalue),inherit.aes = FALSE,size = 3)
wrap_plots(list(p2,p1), ncol=2)
ggsave(paste0(fig.dir,"/TF_activity_by_state_by_clone_final_TFs_chr4.pdf"), height = 5, width = 11)


