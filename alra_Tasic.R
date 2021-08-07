# Imputation analysis of Tasic et. al data

source("../ALRA/alra.R")
source("../ALRA/alraSeurat2.R")
  
library(Seurat) # Seruat Version 2
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Rmagic)
library(ranger)
source('convenience.R')

## Load data

download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746%5Fcells%5Fexon%5Fcounts%2Ecsv%2Egz',   "data/GSE115746_cells_exon_counts.csv.gz")
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746%5Fcomplete%5Fmetadata%5F28706%2Dcells%2Ecsv%2Egz',   "data/GSE115746_complete_metadata_28467-cells.txt.gz")

# exp mat
data_exp <- read.table(
  "data/GSE115746_cells_exon_counts.csv.gz", sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
data_exp <- Matrix(as.matrix(data_exp), sparse = T)

dim(data_exp)

# cell metadata
cell_info <- read.table( "data/GSE115746_complete_metadata_28467-cells.txt.gz", sep = ",", header = T, row.names = 1, stringsAsFactors = F, 
)
data_exp <- data_exp[,colnames(data_exp) %in% cell_info$title]


# create Seurat object
data_original_S <- CreateSeuratObject(raw.data = data_exp, min.genes = 1000, min.cells = 5)
data_original_S <- NormalizeData(data_original_S)
data_original_S@meta.data$cluster <- cell_info[data_original_S@cell.names,"cell_cluster"]
data_original_S <- SetAllIdent(data_original_S, id = "cluster")


# filter low quality cells
# data_original_S <- SubsetData(data_original_S, ident.remove = grep("Low Quality"))

data_original_S


## Imputation

fn ='data/tasic_alra.RDS';
if ( !file.exists(fn)) {
        set.seed(3)
        chose.k <- choose_k(as.matrix(t(data_original_S@data)))
	saveRDS(chose.k,'data/tasic_alra_choose_k.RDS')
        data_alra_S <- alraSeurat2(data_original_S,k=chose.k$k)
	saveRDS(data_alra_S,fn)
}else{
        print(sprintf("Loading %s\n", fn));
	chose.k <- readRDS('data/tasic_alra_choose_k.RDS')
        data_alra_S <- readRDS(fn)
}


# magic
fn ='data/tasic_magic.RDS';
if ( !file.exists(fn)) {
        data_magic <- magic(t(data_original_S@data))
        data_magic_S <- data_original_S
        data_magic_S@data <- Matrix(as.matrix(t(data_magic$result)))
	saveRDS(data_magic_S,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        data_magic_S <- readRDS(fn)
}



# read in dca result
fn <- "data/dca_tasic/mean.tsv"
if (!file.exists(fn)) {
    system.time(system(sprintf('dca data/GSE115746_cells_exon_counts.csv.gz data/dca_tasic')))
}

# DCA
data_dca <- read.table( fn, sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F)
data_dca <- Matrix(as.matrix(data_dca[rownames(data_original_S@data),data_original_S@cell.names]))
data_dca_S <- data_original_S
data_dca_S@data <- data_dca
data_dca_S@raw.data <- data_dca
data_dca_S <- NormalizeData(data_dca_S)
rm(data_dca)

#fofo <- data_dca[rownames(data_original_S@data),data_original_S@cell.names]
#data_dca[1:5,1:5]
#sum(  tolower(data_original_S@cell.names) %in% colnames(data_dca))
#colnames(data_dca)[1:5]
#data_original_S@cell.names[1:5]
#str(data_original_S)
#[1:5]

## Get t-SNE

# function
get.tsne <- function(object, n.pcs = 30){
  object <- ScaleData(object)
#  object <- FindVariableGenes(object, do.plot = F)
#  print(length(object@var.genes))
  object@var.genes <- rownames(object@data)
  object <- RunPCA(object, pcs.compute = n.pcs, do.print = F)
  object <- RunTSNE(
    object, dims.use = 1:n.pcs, seed.use = 1, tsne.method = "FIt-SNE", 
    fast_tsne_path = "/home/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne"
  )
  return(object)
}


# cells to use
cluster_remove <- c(
  grep("Low Quality", unique(data_original_S@ident), value = T),
  grep("Doublet", unique(data_original_S@ident), value = T),
  "Batch Grouping VISp L5 PT Chrna6", "Batch Grouping VISp L5 PT Ctxn3", 
  "High Intronic VISp L5 Endou"
)
cells_use <- data_original_S@cell.names[
  !data_original_S@ident %in% cluster_remove
]


# colors to use
n_cluster <- length(unique(data_original_S@ident)) - length(cluster_remove)
# cols_use <- sample(scales::hue_pal()(n_cluster))
set.seed(3)
cols_use <- sample(rainbow(n = n_cluster, s = 0.6, v = 0.9))


# original
data_original_S <- get.tsne(object = data_original_S, n.pcs = chose.k$k)
tsne_original <- data.frame(data_original_S@dr$tsne@cell.embeddings, ident = data_original_S@meta.data$cluster)[cells_use,]


# ALRA
data_alra_S <- get.tsne(object = data_alra_S, n.pcs = chose.k$k)
tsne_alra <- data.frame(data_alra_S@dr$tsne@cell.embeddings, ident = data_alra_S@meta.data$cluster)[cells_use,]


# MAGIC
data_magic_S <- get.tsne(object = data_magic_S, n.pcs = chose.k$k)
tsne_magic <- data.frame(data_magic_S@dr$tsne@cell.embeddings, ident = data_magic_S@meta.data$cluster)[cells_use,]


# DCA
data_dca_S <- get.tsne(object = data_dca_S, n.pcs = chose.k$k)
tsne_dca <- data.frame(data_dca_S@dr$tsne@cell.embeddings, ident = data_dca_S@meta.data$cluster)[cells_use,]



## OOB error

# function to get random forest model
get.rfr <- function(exprs.mat, Y){
  df <- data.frame(t(exprs.mat), Y = as.factor(Y))
  rfr <- ranger(Y ~ ., df)
  return(rfr)
}


# get model for each method
rfr_original <- get.rfr(
  exprs.mat = as.matrix(data_original_S@data[data_original_S@var.genes,cells_use]), 
  Y = data_original_S@meta.data[cells_use,"cluster"]
)

rfr_alra <- get.rfr(
  exprs.mat = as.matrix(data_alra_S@data[data_alra_S@var.genes,cells_use]), 
  Y = data_alra_S@meta.data[cells_use,"cluster"]
)

rfr_magic <- get.rfr(
  exprs.mat = as.matrix(data_magic_S@data[data_magic_S@var.genes,cells_use]), 
  Y = data_magic_S@meta.data[cells_use,"cluster"]
)

rfr_dca <- get.rfr(
  exprs.mat = as.matrix(data_dca_S@data[data_dca_S@var.genes,cells_use]), 
  Y = data_dca_S@meta.data[cells_use,"cluster"]
)


# use PCs to predict labels
rfr_original <- get.rfr(
  exprs.mat = t(as.matrix(data_original_S@dr$pca@cell.embeddings[cells_use,])), 
  Y = data_original_S@meta.data[cells_use,"cluster"]
)

rfr_alra <- get.rfr(
  exprs.mat = t(as.matrix(data_alra_S@dr$pca@cell.embeddings[cells_use,])), 
  Y = data_alra_S@meta.data[cells_use,"cluster"]
)

rfr_magic <- get.rfr(
  exprs.mat = t(as.matrix(data_magic_S@dr$pca@cell.embeddings[cells_use,])), 
  Y = data_magic_S@meta.data[cells_use,"cluster"]
)

rfr_dca <- get.rfr(
  exprs.mat = t(as.matrix(data_dca_S@dr$pca@cell.embeddings[cells_use,])), 
  Y = data_dca_S@meta.data[cells_use,"cluster"]
)



## Plot

toplot <- 1:nrow(tsne_original)
# add OOB error
g1 <- ggplot(tsne_original[toplot,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, col = ident),shape=16, size = 0.01) + 
  scale_color_manual(values = cols_use) + theme_cowplot() + ggtitle("Original") + 
  theme_tsne+
  annotation_compass(sprintf("Error: %.1f%%", rfr_original$prediction.error*100), "SE")
g2 <- ggplot(tsne_alra[toplot,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, col = ident),shape=16, size = 0.01) + 
  scale_color_manual(values = cols_use) + theme_cowplot() + ggtitle("ALRA") + 
  annotation_compass(sprintf("Error: %.1f%%", rfr_alra$prediction.error*100), "SE") + theme_tsne
g3 <- ggplot(tsne_magic[toplot,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, col = ident),shape=16, size = 0.01) + 
  scale_color_manual(values = cols_use) + theme_cowplot() + ggtitle("MAGIC") + 
           annotation_compass(sprintf("Error: %.1f%%", rfr_magic$prediction.error*100), "SE") + theme_tsne
g4 <- ggplot(tsne_dca[toplot,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, col = ident),shape=16, size =0.01) + scale_color_manual(values = cols_use) + theme_cowplot() + ggtitle("DCA") + 
  annotation_compass(sprintf("Error: %.1f%%", rfr_dca$prediction.error*100), "SE") + theme_tsne
g1<-plot_grid(g1,g2,g3,g4,nrow=1)
g1

ggsave(g1,filename="figs/Tasic-novar-tsne.pdf",width=5.8,height=1.45)





##=============================================##
## Test SVD

fn <- "data/tasic_svd.RDS"
if(!file.exists(fn)){
  set.seed(3)
  alra_result <- alra(as.matrix(t(data_original_S@data)), k = chose.k$k)
  svd_result <- Matrix(as.matrix(t(alra_result$A_norm_rank_k)))
  colnames(svd_result) <- data_original_S@cell.names
  data_svd_S <- data_original_S
  data_svd_S@data <- svd_result
  saveRDS(data_svd_S, file = "data/tasic_svd.RDS")
} else {
  data_svd_S <- readRDS("data/tasic_svd.RDS")
}
data_svd_S <- get.tsne(object = data_svd_S, n.pcs = chose.k$k)
tsne_svd <- data.frame(data_svd_S@dr$tsne@cell.embeddings, ident = data_svd_S@meta.data$cluster)[cells_use,]
rfr_svd <- get.rfr(
  exprs.mat = as.matrix(data_svd_S@data[data_svd_S@var.genes,cells_use]), 
  Y = data_svd_S@meta.data[cells_use,"cluster"]
)
rfr_svd <- get.rfr(
  exprs.mat = t(as.matrix(data_svd_S@dr$pca@cell.embeddings[cells_use,])), 
  Y = data_svd_S@meta.data[cells_use,"cluster"]
)

g1 <- ggplot(tsne_svd[toplot,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, col = ident),shape=16, size = 0.01) + 
  scale_color_manual(values = cols_use) + theme_cowplot() + ggtitle("SVD") + 
  annotation_compass(sprintf("Error: %.1f%%", rfr_svd$prediction.error*100), "SE") + theme_tsne
ggsave(g1,filename="figs/Tasic-novar-tsne-svd.pdf",width=1.45,height=1.45)




