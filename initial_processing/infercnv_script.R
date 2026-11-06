#Anne-Marie Lyne 06/01/2021

#script to run infercnv on Neftel et al (2019) 10x dataset

#wd should be where cell_annotations_file.txt, gene_order_file.txt and 
#gbm_ini_proc.rds are stored

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
out.dir = as.character(args[1]) 
samp = as.character(args[2]) 
pbs_arrayid = as.numeric(args[3]) 
ref = as.character(args[4]) 
samp.file = as.character(args[5])
hmm = as.logical(args[6]) #F
hmm.type = as.character(args[7]) #NA
script.dir = as.character(args[8])
data.dir = as.character(args[9])

source(paste0(script.dir,"/make_cell_type_palette.R"))

clust.alg = 'leiden'
sr.ref = "imm"

setwd(out.dir)

library(Seurat)
library(infercnv)
library(Signac)
options(Seurat.object.assay.version = "v3")

#get list of all samples
sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
                         stringsAsFactors = F, fill = T)
samples = sample.info[,1]
samples.newnames = sample.info[,2]

#read in sample, do preliminary processing on RNA, filter cells, find non-malignant cluster
oligo = readRDS(paste0("../oligo_filt_", samples.newnames[pbs_arrayid], ".rds"))
ct.col = which(colnames(oligo@meta.data)==paste0('cell.type.sr.cl.',sr.ref))

#ensure cell types are combined
oligo@meta.data[,ct.col] = combine_cell_types_imm(oligo@meta.data[,ct.col])


#read in infercnv reference of oligo and mac cells
oligo.nonmalig = readRDS(paste0("../../oligo_nonmalig_merge_",ref,".rds"))


#find all cell types in this sample
all.cell.types = unique(oligo@meta.data[,ct.col])
#find which are potentially malignant (i.e. not oligo or mac)
types.malig = all.cell.types[!(all.cell.types %in% c("Macrophage","Oligo"))]
print(types.malig)

#get cells names for these types
clls = names(oligo$orig.ident)[sapply(oligo@meta.data[,ct.col], function(x) {x %in% types.malig})]
print(length(clls))
oligo = subset(oligo, cells = clls)



#check if any cell types in sample have fewer than two cells and remove if so
Ntype = table(oligo@meta.data[,ct.col])
if(any(Ntype<2)){
  inds = which(Ntype<2)
  for(j in 1:length(inds)){
    inds2 = which(oligo@meta.data[,ct.col]==names(Ntype)[inds[j]])
    oligo = oligo[,-inds2]
  }
}

#check if any cell types in reference have fewer than two cells and remove if so
ct.col.nonmalig = which(colnames(oligo.nonmalig@meta.data)==paste0('cell.type.sr.cl.',sr.ref))
Ntype = table(oligo.nonmalig@meta.data[,ct.col.nonmalig])
if(any(Ntype<2)){
  inds = which(Ntype<2)
  for(j in 1:length(inds)){
    inds2 = which(oligo.nonmalig@meta.data[,ct.col.nonmalig]==names(Ntype)[inds[j]])
    oligo.nonmalig = oligo.nonmalig[,-inds2]
  }
}


#extract counts from both seurat objects
DefaultAssay(oligo) <- "RNA"
DefaultAssay(oligo.nonmalig) <- "RNA"

counts_matrix = GetAssayData(oligo, layer="counts")
counts_matrix_nonmalig = GetAssayData(oligo.nonmalig, layer="counts")

#update cell names with sample
colnames(counts_matrix) = paste(colnames(counts_matrix), oligo$orig.ident, sep = "_")

#merge
counts_matrix2 = merge(counts_matrix, counts_matrix_nonmalig, by="row.names", all = T)
rownames(counts_matrix2) = counts_matrix2[,1]
counts_matrix2 = counts_matrix2[,2:ncol(counts_matrix2)]
counts_matrix2[is.na(counts_matrix2)] = 0


#make annotations_file
cell.types = as.character(c(oligo@meta.data[,ct.col], oligo.nonmalig@meta.data[,ct.col.nonmalig]))
to_save = cbind(colnames(counts_matrix2), cell.types)
write.table(to_save, file = "cell_annotations_file.txt",
            quote = F, sep = "\t", col.names = F, row.names = F)

print(table(cell.types))

#find which types of non-malignant cells are present
non.malig = unique(cell.types)[sapply(unique(cell.types), function(x) {!(x %in% types.malig)})]
print(non.malig)


#run infercnv
options(scipen = 100)
infercnv_obj <- CreateInfercnvObject(counts_matrix2,
                                    paste0(data.dir,"/gene_order_file.txt"),
                                    "cell_annotations_file.txt",
                                    non.malig,
                                    delim = "\t",
                                    chr_exclude = c("chrX", "chrY", "chrM"))

if(hmm){

  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, cutoff=0.1 works well for 10x Genomics
    out_dir=out.dir,
    cluster_by_groups=FALSE,
    cluster_references=FALSE,
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=hmm,
    HMM_type=hmm.type,
    no_prelim_plot=TRUE,
    k_obs_groups = 1,
    analysis_mode='subclusters',
    tumor_subcluster_partition_method=clust.alg
  )

}else{

  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, cutoff=0.1 works well for 10x Genomics
    out_dir=out.dir,
    cluster_by_groups=TRUE,
    cluster_references=FALSE,
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=hmm,
    no_prelim_plot=TRUE,
    no_plot = FALSE,
    k_obs_groups = 1,
    tumor_subcluster_partition_method = "leiden",
    leiden_resolution = 0.01,
    write_expr_matrix = TRUE
  )

}





