## https://www.embopress.org/doi/full/10.15252/embj.20105114

library(Seurat) #V3
library(ggplot2)
library(cowplot)

source("../ALRA/alra.R")
source("convenience.R")


### Process data
# download data from
# https://figshare.com/articles/Single-cell_RNA-Seq_of_human_primary_lung_and_bronchial_epithelium_cells/11981034/1
# and save to data/ folder

# read data
Counts_HBECs <- read.csv('data/Lukassen/Counts_HBECs.csv', stringsAsFactors =F)
Metadata_HBECs <- read.csv('data/Lukassen/Metadata_HBECs.csv', stringsAsFactor = F)


## Seurat
# add metadata
data_S <- CreateSeuratObject(counts = Counts_HBECs, project = "HBEC")
table(data_S@meta.data$orig.ident)
data_S@meta.data$celltype <- Metadata_HBECs[colnames(data_S), "Cell.type"]

# integration
sample_names <- sort(levels(data_S@meta.data$orig.ident))
data_S_list <- lapply(sample_names, function(x){
  subset(data_S, idents = x)
})

data_anchors <- FindIntegrationAnchors(object.list = data_S_list)
data_S <- IntegrateData(anchorset = data_anchors)
# data_S <- NormalizeData(data_S)
data_S <- ScaleData(data_S)
# data_S <- FindVariableFeatures(data_S)
data_S <- RunPCA(data_S)
data_S <- RunUMAP(data_S, dims = 1:30)
UMAPPlot(data_S, group.by = "orig.ident")
UMAPPlot(data_S, group.by = "celltype")



## ALRA
DefaultAssay(data_S) <- "RNA"
set.seed(3)
data_S <- RunALRA(data_S, setDefaultAssay = F)



## Label ACE2+/FURIN+/TMPRSS2+ cells
data_S@meta.data$triple.pos <- as.numeric(colSums(data_S@assays$RNA@data[c("ACE2","FURIN","TMPRSS2"),] > 0) == 3)
table(data_S@meta.data$triple.pos)

data_S@meta.data$triple.pos.alra <- as.numeric(colSums(data_S@assays$alra@data[c("ACE2","FURIN","TMPRSS2"),] > 0) == 3)
table(data_S@meta.data$triple.pos.alra)



## Plot

gg1 <- ggplot(data = data.frame(
  Dim1 = data_S@reductions$umap@cell.embeddings[,1], 
  Dim2 = data_S@reductions$umap@cell.embeddings[,2],
  celltype = data_S@meta.data$celltype
)) + geom_point(aes(Dim1, Dim2, col = celltype), shape = 16, size = 0.1) + ggtitle("A)") + 
  theme_tsne + theme(plot.title = element_text(face = "bold"))
ggsave(gg1, filename = "figs/lukassen_celltype.pdf", width = 1.9, height = 1.95)
ggsave(get_legend(gg1 + guides(color = guide_legend(override.aes = list(size = 3))) + 
                    theme(legend.position="right", legend.text=element_text(size=8), legend.key.height=unit(0.25,"cm"), legend.title=element_blank())), 
       filename = "figs/lukassen_celltype_legend.pdf", width = 1.2, height = 1.9)

gg2 <- ggplot(data = data.frame(
  Dim1 = data_S@reductions$umap@cell.embeddings[,1], 
  Dim2 = data_S@reductions$umap@cell.embeddings[,2]
)) + geom_point(aes(Dim1, Dim2), size = 0.1, col = "#CFCFCF") + 
  geom_point(data = data.frame(
    Dim1 = data_S@reductions$umap@cell.embeddings[data_S@meta.data$triple.pos==1,1], 
    Dim2 = data_S@reductions$umap@cell.embeddings[data_S@meta.data$triple.pos==1,2]
  ), aes(Dim1,Dim2), col = "#A81E22", size = 0.25) + ggtitle("B)") + theme_tsne + theme(plot.title = element_text(face = "bold"))
ggsave(gg2, filename = "figs/lukassen_original.pdf", width = 1.9, height = 1.95)

gg3 <- ggplot(data = data.frame(
  Dim1 = data_S@reductions$umap@cell.embeddings[,1], 
  Dim2 = data_S@reductions$umap@cell.embeddings[,2]
)) + geom_point(aes(Dim1, Dim2), size = 0.1, col = "#CFCFCF") + 
  geom_point(data = data.frame(
    Dim1 = data_S@reductions$umap@cell.embeddings[data_S@meta.data$triple.pos.alra==1,1], 
    Dim2 = data_S@reductions$umap@cell.embeddings[data_S@meta.data$triple.pos.alra==1,2]
  ), aes(Dim1,Dim2), col = "#A81E22", size = 0.25) + ggtitle("C)") + theme_tsne + theme(plot.title = element_text(face = "bold"))
ggsave(gg3, filename = "figs/lukassen_alra.pdf", width = 1.9, height = 1.95)



## Test different scale factors

data_S_factor <- list()
for(sf in c(5000,20000,10000,1000)){
  data_S_factor[[as.character(sf)]] <- NormalizeData(data_S, scale.factor = sf)
  data_S_factor[[as.character(sf)]]@tools$RunALRA <- NULL
  data_S_factor[[as.character(sf)]]@assays$alra <- NULL
  set.seed(3)
  data_S_factor[[as.character(sf)]] <- RunALRA(data_S_factor[[as.character(sf)]], setDefaultAssay = F, q.k = 10)
  data_S_factor[[as.character(sf)]]@meta.data$triple.pos <- 
    as.numeric(colSums(data_S_factor[[as.character(sf)]]@assays$RNA@data[c("ACE2","FURIN","TMPRSS2"),] > 0) == 3)
  print(table(data_S_factor[[as.character(sf)]]@meta.data$triple.pos))
  
  data_S_factor[[as.character(sf)]]@meta.data$triple.pos.alra <- 
    as.numeric(colSums(data_S_factor[[as.character(sf)]]@assays$alra@data[c("ACE2","FURIN","TMPRSS2"),] > 0) == 3)
  print(table(data_S_factor[[as.character(sf)]]@meta.data$triple.pos.alra))
}



