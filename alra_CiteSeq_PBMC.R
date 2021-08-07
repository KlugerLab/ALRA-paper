##=======================================================##
## CITE-Seq PBMC
## ALRA
##=======================================================##


# load Seurat v3.0
devtools::load_all("../seurat/")
set.seed(3)

library(Matrix)
library(ggplot2)
library(reshape2)
library(cowplot)



## Prepare data

# read in data
# get ADT data
fn <-  "./data/GSE100866_PBMC_vs_flow_10X-RNA_umi.csv.gz"
if (!file.exists(fn) ){
    download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FPBMC%5Fvs%5Fflow%5F10X%2DRNA%5Fumi%2Ecsv%2Egz',fn)
}
data_exp_rna <- read.table(
fn, 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
# filter cells
human_ratio_sum_per_cell <- colSums(data_exp_rna[grep("HUMAN_", rownames(data_exp_rna), fixed = T),]) / colSums(data_exp_rna)
human_cell_idx <- which(human_ratio_sum_per_cell >= 0.9)
data_exp_rna <- data_exp_rna[,human_cell_idx]


# filter mouse genes
data_exp_rna <- data_exp_rna[grep("HUMAN_", rownames(data_exp_rna), fixed = T),]
dim(data_exp_rna)

# filter low-expressed human genes
num_of_cells_in_gene <- rowSums(data_exp_rna > 0)
data_exp_rna <- data_exp_rna[num_of_cells_in_gene > 10,]

# change format
data_exp_rna <- Matrix(as.matrix(data_exp_rna), sparse = T)
mean(data_exp_rna != 0)

# get ADT data
fn <-  "./data/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv.gz"
if (!file.exists(fn) ){
    download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FPBMC%5Fvs%5Fflow%5F10X%2DADT%5Fumi%2Ecsv%2Egz',fn)
}
data_exp_adt <- read.table(
fn, 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
data_exp_adt <- data_exp_adt[,human_cell_idx]

fn <- "./data/GSE100866_PBMC_vs_flow_10X-ADT_clr-transformed.csv.gz"
if (!file.exists(fn) ){
    download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FPBMC%5Fvs%5Fflow%5F10X%2DADT%5Fclr%2Dtransformed%2Ecsv%2Egz',fn)
}

# CLR-tranformed data
data_norm_adt <- read.table(
fn, 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
data_norm_adt <- data_norm_adt[,human_cell_idx]


# Seurat object
data_S <- CreateSeuratObject(counts = data_exp_rna)
data_S[["ADT"]] <- CreateAssayObject(counts = as.matrix(data_norm_adt))

data_S <- NormalizeData(data_S)
#data_S <- NormalizeData(data_S, assay = "ADT", normalization.method = "CLR")
data_S@assays$ADT@data <- as.matrix(data_norm_adt)





## Sort cells

# k-means
kmeans_test <- kmeans(x = GetAssayData(data_S,assay="ADT")["CD3",], centers = 2)
data_S$kmeans_test <- as.factor(kmeans_test$cluster)

# cluster based on CD3,CD19,CD14
genes_to_gate <- c("CD3","CD19","CD14")
for(gene in genes_to_gate){
  kmeans_gene <- kmeans(x = GetAssayData(data_S,assay="ADT")[gene,], centers = 2)
  kmeans_gene <- kmeans_gene$cluster
  # determine which is gene+/gene-
  if(
    max(GetAssayData(data_S,assay="ADT")[gene,kmeans_gene == 1]) > 
    max(GetAssayData(data_S,assay="ADT")[gene,kmeans_gene == 2])
  ){
    kmeans_gene[kmeans_gene == 1] <- paste(gene,"+",sep = "")
    kmeans_gene[kmeans_gene == 2] <- paste(gene,"-",sep = "")
  } else {
    kmeans_gene[kmeans_gene == 2] <- paste(gene,"+",sep = "")
    kmeans_gene[kmeans_gene == 1] <- paste(gene,"-",sep = "")
  }
  data_S@meta.data[,gene] <- as.factor(kmeans_gene)
}

#FeatureScatter(
#  object = data_S, feature1 = paste("adt_",gene1,sep = ""), feature2 = paste("adt_",gene2,sep = ""),
#  group.by = "CD3"
#)

#pdf("./gating_results.pdf")
#for(gene1 in genes_to_gate[1:2]){
#  for(gene2 in genes_to_gate[-1]){
#    if(gene1 == gene2){next}
#    data_S@meta.data$double <- paste(data_S@meta.data[,gene1], data_S@meta.data[,gene2], sep = ".")
#    print(FeatureScatter(
#      object = data_S, feature1 = paste("adt_",gene1,sep = ""), feature2 = paste("adt_",gene2,sep = ""),
#      group.by = "double"
#    ))
#  }
#}
#dev.off()


# add cell type
data_S@meta.data$celltype <- "Other"
data_S@meta.data$celltype[
  data_S@meta.data$CD3 == "CD3+" & data_S@meta.data$CD19 == "CD19-" & data_S@meta.data$CD14 == "CD14-"
] <- "T_cells"
data_S@meta.data$celltype[
  data_S@meta.data$CD3 == "CD3-" & data_S@meta.data$CD19 == "CD19+" & data_S@meta.data$CD14 == "CD14-"
] <- "B_cells"
data_S@meta.data$celltype[
  data_S@meta.data$CD3 == "CD3-" & data_S@meta.data$CD19 == "CD19-" & data_S@meta.data$CD14 == "CD14+"
] <- "Monocytes"
table(data_S@meta.data$celltype)

#pdf("./gating_labels.pdf")
#for(gene1 in genes_to_gate[1:2]){
#  for(gene2 in genes_to_gate[-1]){
#    if(gene1 == gene2){next}
#    print(FeatureScatter(
#      object = data_S, feature1 = paste("adt_",gene1,sep = ""), feature2 = paste("adt_",gene2,sep = ""),
#      group.by = "celltype"
#    ))
#  }
#}
#dev.off()


# gating markers on RNA level
genes_to_gate_rna <- c("CD3"="CD3D","CD19"="CD19","CD14"="CD14")
genes_to_gate_rna <- paste("HUMAN-", genes_to_gate_rna, sep = "")
names(genes_to_gate_rna) <- genes_to_gate







## Bulk RNA-seq data

# read in data
PBMC_bulk <- read.table(
  "./data/GSE64655_Normalized_transcript_expression_in_human_immune_cells.txt.gz",
  sep = "\t", header = T, stringsAsFactors = F, skip = 4, quote="\""
)
PBMC_bulk$Gene.Symbol.H <- paste("HUMAN-", gsub("_","-",PBMC_bulk$Gene.Symbol), sep = "")


# filter bulk data
PBMC_bulk <- subset(PBMC_bulk, Gene.Symbol.H %in% rownames(data_S))


# match bulk label & scRNA label, in SAME order
bulklabel2do <- list(c("HD30_B_d0","HD31_B_d0"),
                     c("HD30_Mono_d0","HD31_Mono_d0"),
                     c("HD30_T_d0","HD31_T_d0"))
sclabel2do <- c("B_cells", "Monocytes", "T_cells")
n_label2do <- length(sclabel2do)


# get zero genes for each label
getZeroGenes <- function(bulk.data, sc.data, bulk.genes = NULL){
  if(is.null(bulk.genes)){
    bulk.genes <- rownames(bulk.data)
  }
  zero.genes <- unique(bulk.genes[rowSums(bulk.data) == 0])
  n.zero.genes <- length(zero.genes)
  # deal with gene symbols with multiple IDs
  zero.idx.rm <- integer()
  for(ii in 1:n.zero.genes){
    zero.gene <- zero.genes[ii]
    if(sum(bulk.data[bulk.genes %in% zero.gene,]) > 0){
      zero.idx.rm <- c(zero.idx.rm, ii)
      cat("Remove gene ", zero.gene, "\n", sep = "")
    }
  }
  if(length(zero.idx.rm) > 0){
    zero.genes <- zero.genes[-zero.idx.rm]
  } 
  cat("# zero genes: ", length(zero.genes), "\n", sep = "")
  cat("# zero entries: ", sum(sc.data[zero.genes,] == 0), "\n", sep = "") 
  return(zero.genes)
}

PBMC_bulk_zero_genes <- list()
for(i in 1:n_label2do){
  PBMC_bulk_zero_genes[[i]] <- getZeroGenes(
    bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
    sc.data = GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]], 
    bulk.genes = PBMC_bulk$Gene.Symbol.H
  )
}
names(PBMC_bulk_zero_genes) <- sclabel2do


# double check zero genes
for(i in 1:n_label2do){
  print(sum(PBMC_bulk[PBMC_bulk$Gene.Symbol.H %in% PBMC_bulk_zero_genes[[i]],bulklabel2do[[i]]]))
}

which(sapply(PBMC_bulk_zero_genes[[i]], FUN = function(x, dat, genes){
  sum(dat[genes %in% x,])
}, dat = PBMC_bulk[,bulklabel2do[[i]]], genes = PBMC_bulk$Gene.Symbol.H) > 0)


# original zero ratio for zero genes
for(i in 1:n_label2do){
  cat(
    mean(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0), 
    "\n", sep = ""
  )
}


# write zero genes in table
write.table(
  PBMC_bulk[PBMC_bulk$Gene.Symbol.H %in% unique(unlist(PBMC_bulk_zero_genes)), 
            c("Gene.ID", "Gene.Symbol", unlist(bulklabel2do))],
  file = "data/PBMC_bulk_day0_zero_genes", sep = "\t", row.names = F, col.names = T, quote = F)





## Imputation

# ALRA
set.seed(3)
data_S <- RunALRA(data_S)
options(error=traceback)


# MAGIC
library(Rmagic)
magic_result <- magic(t(GetAssayData(data_S, assay = "RNA", slot = "data")))
save(list = "magic_result", file = "data/cite_pbmc_magic.RData")

magic_norm <- t(magic_result$result)


# SAVER
library(SAVER)
fn <- 'data/cite_pbmc_saver.RData'
if ( !file.exists(fn)) {
        print(sprintf("Generating %s\n", fn));
        saver_result <- saver( GetAssayData(data_S, assay = "RNA", slot = "counts"))
        save(list = c("saver_result"), file =fn)
}else{
    print(sprintf("Loading %s\n", fn));
    load(fn)
}
saver_norm <- sample.saver(saver_result, seed = 3)


# scImpute

#scimpute_norm <- Matrix(as.matrix(scimpute_norm), sparse = T)

# scImpute
fn <- 'data/cite_pbmc_scimpute.RData'
if ( !file.exists(fn)) {
        print(sprintf("Generating %s\n", fn));
        library(scImpute)
        write.table( data_exp_rna, "./GSE100866_PBMC_rna.csv", sep = ",", quote = F, row.names = T, col.names = T)
#        scimpute_result <- scimpute( count_path = "GSE100866_PBMC_rna.csv",
#          out_dir = "./data/scimpute_cite_pbmc_0.01/", 
#          Kcluster = 5, 
#          drop_thre = 0.01)
#        scimpute_norm0.01 <- read.table( "./data/scimpute_cite_pbmc_0.01/scimpute_count.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
#
        scimpute_result <- scimpute( count_path = "GSE100866_PBMC_rna.csv",
          out_dir = "./data/scimpute_cite_pbmc/", 
          Kcluster = 5)
        scimpute_norm_default <- read.table( "data/scimpute_cite_pbmc/scimpute_count.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
        rownames(scimpute_norm_default) <- gsub("_","-", rownames(scimpute_norm_default))
        save(list = c("scimpute_norm_default"), file =fn)
}else{
    print(sprintf("Loading %s\n", fn));
    load(fn)
}



## Evaluation

# overall completion performance
mean(GetAssayData(data_S, assay = "RNA", slot = "counts") == 0)
mean(GetAssayData(data_S, assay = "alra", slot = "data") == 0)
mean(magic_norm <= 0)
mean(saver_norm == 0)
mean(scimpute_norm_default == 0)
#mean(dca_norm <= 0)


# zero preservation
zero_preserve_alra <- rep(0, n_label2do)
names(zero_preserve_alra) <- sclabel2do
zero_preserve_magic <- zero_preserve_alra
zero_preserve_saver <- zero_preserve_alra
zero_preserve_scimpute <- zero_preserve_alra
for(i in 1:n_label2do){
  zero_preserve_alra[i] <-
    sum(GetAssayData(data_S, assay = "alra")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) /
    sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0)
  zero_preserve_magic[i] <- sum(magic_norm[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] <= 0) /
    sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) 
  zero_preserve_saver[i] <- sum(saver_norm[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) /
    sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0)
  zero_preserve_scimpute[i] <- sum(scimpute_norm_default[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) / 
    sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0)
}
print(zero_preserve_scimpute)


# (zero preservation ratio - sparsity) per cell type
zero_complete_alra <- rep(0, n_label2do)
names(zero_complete_alra) <- sclabel2do
zero_complete_magic <- zero_complete_alra
zero_complete_saver <- zero_complete_alra
zero_complete_scimpute <- zero_complete_alra
for(i in 1:n_label2do){
   zero_complete_alra[i] <- mean(
    (GetAssayData(data_S, assay = "alra")[,data_S$celltype == sclabel2do[i]] > 0) & 
    (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0))
   zero_complete_magic[i] <- mean(
    magic_norm[,data_S$celltype == sclabel2do[i]] > 0 & 
    (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0))
   zero_complete_saver[i] <- mean(saver_norm[,data_S$celltype == sclabel2do[i]] > 0 &
   (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0))
  zero_complete_scimpute[i] <-  mean( scimpute_norm_default[,data_S$celltype == sclabel2do[i]] > 0 & (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0)
  )
}

for(i in 1:3){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_alra[i], zero_complete_alra[i]))
}

for(i in 1:3){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_magic[i], zero_complete_magic[i]))
}

for(i in 1:3){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_saver[i], zero_complete_saver[i]))
}

for(i in 1:3){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_scimpute[i], zero_complete_scimpute[i]))
}


# Check sensitivity to the choice of $k$
K_s <- 15:25
for (K_i in 1:length(K_s)  ) {
    K_ <- K_s[K_i]
    DefaultAssay(data_S) <- "RNA"
    data_S <- RunALRA(data_S,k=K_)
    for(i in 1:n_label2do){
        zero_preserve <- sum(GetAssayData(data_S, assay = "alra")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) / sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0)
       zero_complete <- mean( (GetAssayData(data_S, assay = "alra")[,data_S$celltype == sclabel2do[i]] > 0) & (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0))
       print(sprintf('%d, %s: %.2f, %.2f', K_, sclabel2do[i], zero_preserve, zero_complete))
    }
}





################################
# Test scale factor
################################

for(sf in c(1000,5000,10000,20000)){
  # normalize with different size factor
  data_S <- NormalizeData(data_S, scale.factor = sf)
  
  # ALRA
  set.seed(3)
  data_S <- RunALRA(data_S, setDefaultAssay = F)
  zero_preserve_alra <- rep(0, n_label2do)
  names(zero_preserve_alra) <- sclabel2do
  for(i in 1:n_label2do){
    zero_preserve_alra[i] <-
      sum(GetAssayData(data_S, assay = "alra")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0) /
      sum(GetAssayData(data_S, assay = "RNA")[PBMC_bulk_zero_genes[[i]],data_S$celltype == sclabel2do[i]] == 0)
  }
  zero_complete_alra <- rep(0, n_label2do)
  names(zero_complete_alra) <- sclabel2do
  for(i in 1:n_label2do){
    zero_complete_alra[i] <- mean(
      (GetAssayData(data_S, assay = "alra")[,data_S$celltype == sclabel2do[i]] > 0) & 
        (GetAssayData(data_S, assay = "RNA")[,data_S$celltype == sclabel2do[i]] == 0))
  }
  for(i in 1:3){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_alra[i], zero_complete_alra[i]))
  }
  data_S@assays$alra <- NULL
  data_S@tools$RunALRA <- NULL
}


