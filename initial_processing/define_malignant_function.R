define_malignant_cells <- function(oligo, cl.res, icnv.diff, sr.ref, pq.loc, iCNV.ref, iCNV.obs, data.dir){
  
  
  #define SingleR cell types which are generally not malignant
  cell.types.nonmalig = c("Oligo","Macrophage","GABAergic","Glutamatergic","T-cells","B-cells")
  
  #find which meta data column with SingleR cell type
  ct.col = which(colnames(oligo@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
  oligo@meta.data[,ct.col] = combine_cell_types_imm(oligo@meta.data[,ct.col])
  
  #find if astro or oligo
  cancer = unlist(strsplit(oligo$orig.ident[1],"_",fixed=T))[1]
  md.col = which(colnames(oligo@meta.data)==paste0("SCT_snn_res.",cl.res))
  cl.uni = unique(oligo@meta.data[,md.col])
  if(cancer=="oligo"){
    #relatively simple for IDH-O
    #see if 1p deletion present
    malig.1p7q.stat = sapply(cl.uni, function(x) {
      
      #if in cell.types.nonmalig, call as non-malignant
      #otherwise, compute 1p median and ddfine as malignant if <0.975
      cell.type = unique(oligo@meta.data[which(oligo@meta.data[,md.col]==x),ct.col])
      if(cell.type %in% cell.types.nonmalig){
        return(FALSE)
      }else{
        #compute median of 1p infercnv value
        med.1p = median(oligo$p1.av[which(oligo@meta.data[,md.col]==x)], na.rm = T)
        if(med.1p<(1-icnv.diff)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    })
  }else{
    #more complicated for IDH-A where we don't have a specific CNA present
    #some samples have clear (non 7q) CNA
    if(oligo$orig.ident[1] %in% c("astro_P0027_rec1_S0053","astro_P0027_ini_S0052",
                                  "astro_P0015_ini_S0028","astro_P0028_ini_S0054",
                                  "astro_P0028_rec1_S0055","astro_P0018_rec1_S0035",
                                  "astro_P0019_rec1_S0037")){
      if(oligo$orig.ident[1] %in% c("astro_P0027_rec1_S0053","astro_P0027_ini_S0052","astro_P0015_ini_S0028")){
        #use chr 21
        #get all genes for alternative chromosome
        altchr.genes = gene.order[gene.order[,2]=="chr21",1]
        #average inferCNV for chr21 genes in reference cells
        altchr.ref.av = colMeans(iCNV.ref[rownames(iCNV.ref) %in% altchr.genes,])
        #average for observation cells
        altchr.obs.av = colMeans(iCNV.obs[rownames(iCNV.obs) %in% altchr.genes,])
      }else if(oligo$orig.ident[1] %in% c("astro_P0028_ini_S0054","astro_P0028_rec1_S0055")){
        #use chr 12q
        altchr.rows = which(pq.loc[,1]=="chr12" & grepl("^q",pq.loc[,4]))
        altchr.st = pq.loc[altchr.rows[1],2] #start of 7q
        altchr.en = pq.loc[altchr.rows[length(altchr.rows)],3] #end of 7q
        #get all genes in this region on 7q
        altchr.genes = gene.order[gene.order[,2]=="chr12" & gene.order[,3]>altchr.st & gene.order[,4]<altchr.en,1]
        #average inferCNV for 7q genes in reference cells
        logi = rownames(iCNV.ref) %in% altchr.genes
        #find number of 7q genes
        naltchr = sum(logi)
        indsaltchr = which(logi)
        altchr.ref.av = colMeans(iCNV.ref[logi,])
        #average for observation cells
        altchr.obs.av = colMeans(iCNV.obs[logi,])
      }else if(oligo$orig.ident[1]=="astro_P0018_rec1_S0035"){
        #use chr 15
        #get all genes
        altchr.genes = gene.order[gene.order[,2]=="chr15",1]
        #average inferCNV for chr15 genes in reference cells
        altchr.ref.av = colMeans(iCNV.ref[rownames(iCNV.ref) %in% altchr.genes,])
        #average for observation cells
        altchr.obs.av = colMeans(iCNV.obs[rownames(iCNV.obs) %in% altchr.genes,])
      }else if(oligo$orig.ident[1]=="astro_P0019_rec1_S0037"){
        #use chr 4q
        altchr.rows = which(pq.loc[,1]=="chr4" & grepl("^q",pq.loc[,4]))
        altchr.st = pq.loc[altchr.rows[1],2] #start of 4q
        altchr.en = pq.loc[altchr.rows[length(altchr.rows)],3] #end of 4q
        #get all genes in this region on 4q
        altchr.genes = gene.order[gene.order[,2]=="chr4" & gene.order[,3]>altchr.st & gene.order[,4]<altchr.en,1]
        #average inferCNV for 4q genes in reference cells
        logi = rownames(iCNV.ref) %in% altchr.genes
        #find number of 4q genes
        naltchr = sum(logi)
        indsaltchr = which(logi)
        altchr.ref.av = colMeans(iCNV.ref[logi,])
        #average for observation cells
        altchr.obs.av = colMeans(iCNV.obs[logi,])
      }
      altchr.av = sapply(paste(colnames(oligo@assays$SCT),
                               oligo$orig.ident[1], sep = "_"),
                         function(y) {
                           if(y %in% names(c(altchr.obs.av,altchr.ref.av))){
                             c(altchr.obs.av,altchr.ref.av)[which(names(c(altchr.obs.av,altchr.ref.av))==y)]
                           }else{
                             NA
                           }
                         })
      names(altchr.av) = colnames(oligo@assays$SCT)
      oligo$altchr.av = altchr.av
      #boxplot showing alternative chromsome per cluster
      df = data.frame(cluster = oligo@meta.data[,md.col], altchr = oligo$altchr.av,
                      cell.type = oligo@meta.data[,ct.col])
      p5 = ggplot(df, aes(x=cluster, y=altchr, fill=cell.type)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 0.05, position=position_jitter(width=0.2,height=0)) +
        geom_hline(yintercept=(1-icnv.diff), color = "blue", linewidth=1.5) +
        geom_hline(yintercept=(1+icnv.diff), color = "blue", linewidth=1.5)
      print(p5)
      ggsave(paste0(data.dir, "/figures/boxplot_altchr_",oligo$orig.ident[1],".png"), width=18, height=13)
      
      malig.1p7q.stat = sapply(cl.uni, function(x) {
        
        #first check cell type for this cluster
        #if in non-malig list, call as non-malignant
        #otherwise, compute altchr median and define as malignant if <thresh
        cell.type = unique(oligo@meta.data[which(oligo@meta.data[,md.col]==x),ct.col])
        if(cell.type %in% cell.types.nonmalig){
          return(FALSE)
        }else{
          #compute median of altchr infercnv value
          med.altchr = median(oligo$altchr.av[which(oligo@meta.data[,md.col]==x)],na.rm=T)
          if(med.altchr>(1+icnv.diff) || med.altchr<(1-icnv.diff)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      })
    }else{
      #for remaining IDH-A, check 7q amplification
      malig.1p7q.stat = sapply(cl.uni, function(x) {
        
        #check if cluster is in non-malig list, call as non-malignant
        #otherwise, compute 7q median and define as malignant if <0.975
        cell.type = unique(oligo@meta.data[which(oligo@meta.data[,md.col]==x),ct.col])
        if(cell.type %in% cell.types.nonmalig){
          return(FALSE)
        }else{
          #compute median of 1p infercnv value
          med.7q = median(oligo$q7.av[which(oligo@meta.data[,md.col]==x)], na.rm=T)
          if(med.7q>(1+icnv.diff)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      })
    }
  }
  
  #not all samples have CNA, so also check if clusters have SNV
  malig.SNV.stat = sapply(cl.uni, function(x) {
    cell.type = unique(oligo@meta.data[which(oligo@meta.data[,md.col]==x),ct.col])
    if(cell.type %in% cell.types.nonmalig){
      return(FALSE)
    }else{
      SNV.rel = oligo$IDH.mut[oligo@meta.data[,md.col]==x]
      sum(!(SNV.rel %in% c("non_mutant","not_expressed")))>0
    }
  })
  
  #for each cell, find which cluster and define as malignant if either 1p/7q or SNV present (astro)
  #just based on 1p if oligo
  iCNV.malig.ident = sapply(seq(ncol(oligo)), function(x) {
    clust = oligo@meta.data[x,md.col]
    if(grepl("oligo",oligo$orig.ident[1])){
      if(malig.1p7q.stat[which(cl.uni==clust)]){
        return("malignant")
      }else{
        return("non-malig")
      }
    }else{
      if(malig.1p7q.stat[which(cl.uni==clust)] | malig.SNV.stat[which(cl.uni==clust)]){
        return("malignant")
      }else{
        return("non-malig")
      }
    }
  })
  names(iCNV.malig.ident) = names(oligo$orig.ident)
  oligo$iCNV.malig.ident = iCNV.malig.ident
  
  #works reasonably well for most samples, but need some changes for specific populations in some samples
  if(oligo$orig.ident[1]=="astro_P0015_rec1_S0029"){
    #some glutamatergic neurons which cluster with malignant cells on UMAP have CNV
    emb.umap = Embeddings(oligo, reduction = "umap.rna")
    oligo$iCNV.malig.ident[which(oligo$cell.type.sr.cl.imm=="Glutamatergic" &
                                   emb.umap[,1]>(-3) & emb.umap[,1]<8 &
                                   emb.umap[,2]>(-7) & emb.umap[,2]<0)] = "malignant"
  }
  if(oligo$orig.ident[1]=="astro_P0018_ini_S0034"){
    #add cluster
    #hard to tell since no clear CNA, but OPC and astro clusters have IDH1 mutation so define as malignant
    clls = names(oligo$orig.ident)[which(oligo$cell.type.sr.cl.imm=="OPC" | oligo$cell.type.sr.cl.imm=="Astro")]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="astro_P0019_ini_S0036"){
    #add clusters
    #all OPC and astro malignant by last 50 on chr19q in ATAC
    #not sure about oligo cluster which is attached to astro cluster on umap, leave as non-malignant
    clls = names(oligo$orig.ident)[which(oligo$cell.type.sr.cl.imm=="OPC" | oligo$cell.type.sr.cl.imm=="Astro")]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="astro_P0019_rec1_S0037"){
    #add cluster
    #annotated GABAergic, but clearly malignant by CNA
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==7)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="astro_P0029_ini_S0056"){
    #add clusters
    #hard to tell, but two of the astro clusters malignant by CNA
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==10 | oligo@meta.data[,md.col]==14)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0002_rec1_S0005"){
    #add cluster
    #annotated B-cells but clearly malignant by 1p/19q
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==10)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0004_rec1_S0010"){
    #add cluster
    #annotated B-cells but clearly malignant by 1p/19q
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==12)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0005_rec1_S0012"){
    #add cluster 13, GABAergic but clearly malignant by 1p
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==12)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0008_rec1_S0018"){
    #add OPC cluster as seems to have IDH1 mutation
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==3)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0012_rec1_S0024"){
    #add clusters 2 and 5 OPC
    #just misses threshold, and seems malignant by ATAC
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==2 |
                                           oligo@meta.data[,md.col]==4)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0022_ini_S0042"){
    #set all Astro, Oligo and OPC as malignant, seems that way when run infercnv atac
    #with just macrophage and endo as reference cells
    clls = names(oligo$orig.ident)[which(oligo$cell.type.sr.cl.imm=="OPC" |
                                           oligo$cell.type.sr.cl.imm=="Astro"|
                                           oligo$cell.type.sr.cl.imm=="Oligo")]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0024_ini_S0046"){
    #add astro cluster as seems malignant by 19q in ATAC
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==13)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  if(oligo$orig.ident[1]=="oligo_P0031_rec1_S0061"){
    #add OPC cluster as seems malignant by 19q in ATAC and IDH mutation
    clls = names(oligo$orig.ident)[which(oligo@meta.data[,md.col]==18)]
    oligo$iCNV.malig.ident[which(names(oligo$iCNV.malig.ident) %in% clls)] = "malignant"
  }
  return(oligo)
  
}






