#source("https://raw.githubusercontent.com/KlugerLab/ALRA/quantiles/alra.R")
source("../ALRA/alra.R")

library(Seurat) # Version 2.3.4
library(Rmagic)
library(SAVER)
library(scImpute)
library(gridExtra)



## Data preprocessing

# function to read 10x data
read10X <- function(folder = NULL, matFile = NULL, geneFile = NULL, cellFile = NULL, 
                    suffix = "", sep = "_", gz = T){
  if(!is.null(folder)){
    if(gz){
      matFile <- paste(folder, "/matrix.mtx.gz", sep = "")
      geneFile <- paste(folder, "/genes.tsv.gz", sep = "")
      cellFile <- paste(folder, "/barcodes.tsv.gz", sep = "")
    } else {
      matFile <- paste(folder, "/matrix.mtx", sep = "")
      geneFile <- paste(folder, "/genes.tsv", sep = "")
      cellFile <- paste(folder, "/barcodes.tsv", sep = "")
    }
  }
  
  geneNames <- read.table(geneFile, header = F, sep = "\t", as.is = T)[,2]
  cellNames <- paste(read.table(cellFile, header = F, sep = "\t", as.is = T)[,1], suffix, sep = sep)
  
  # add suffix to duplicate gene names
  if(max(table(geneNames)) > 1){
    for(dupgene in names(which(table(geneNames) != 1))){
      geneidx <- which(geneNames == dupgene)
      for(ii in 2:length(geneidx)){
        geneNames[geneidx[ii]] <- paste(dupgene, ii-1, sep = ".")
      }
    }
  }
  
  rawMat <- readMM(matFile)
  rownames(rawMat) <- geneNames
  colnames(rawMat) <- cellNames
  
  return(rawMat)
}


# download data
if ( !file.exists('data/derm_S.rds')){
  # download data
  urls.download <- c(
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122043/suppl/GSE122043_genes.tsv.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453537/suppl/GSM3453537_e14Control_barcodes.tsv.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453537/suppl/GSM3453537_e14Control_matrix.mtx.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453538/suppl/GSM3453538_e14Control_replicate_barcodes.tsv.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453538/suppl/GSM3453538_e14Control_replicate_matrix.mtx.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453540/suppl/GSM3453540_e14LOF_barcodes.tsv.gz",
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3453nnn/GSM3453540/suppl/GSM3453540_e14LOF_matrix.mtx.gz"
  )
  names(urls.download) <- sapply(urls.download, FUN = function(x){
    x.str <- unlist(strsplit(x, split = "/", fixed = T))
    x.l <- length(x.str)
    return(x.str[x.l])
  })
  for(i in seq_along(urls.download)){
    download.file(
      url = urls.download[i], 
      destfile = paste("data/",names(urls.download)[i], sep = "")
    )
  }
  
  # load data
  sample_names <- c(
    "GSM3453537_e14Control","GSM3453538_e14Control_replicate","GSM3453540_e14LOF"
  )
  n_sample <- length(sample_names)
  
  data_list <- list()
  
  for(i in 1:n_sample){
    data_list[[i]] <- read10X(
      matFile = paste("./data/", sample_names[i], "_matrix.mtx.gz", sep = ""), 
      cellFile = paste("./data/", sample_names[i], "_barcodes.tsv.gz", sep = ""),
      geneFile = "./data/GSE122043_genes.tsv.gz", 
      suffix = sample_names[i], sep = "-"
    )
  }
  names(data_list) <- sample_names
  
  # Seurat object for each sample
  data_S_list <- lapply(data_list, function(x){
    x_S <- CreateSeuratObject(raw.data = x, min.genes = 1000, names.delim = "-", names.field = 3)
    x_S <- NormalizeData(x_S)
    x_S <- ScaleData(x_S)
    x_S <- FindVariableGenes(x_S, do.plot = F)
    x_S <- RunPCA(
      x_S, pc.genes = rownames(subset(x_S@hvg.info, gene.dispersion > 0.8)),
      pcs.compute = 10, do.print = F
    )
    x_S <- RunTSNE(x_S, dims.use = 1:10)
    x_S <- FindClusters(x_S, dims.use = 1:10, resolution = 0.1, print.output = F)
    
    return(x_S)
  })
  names(data_S_list) <- sample_names
  
  # select only dermal cells based on Col1a1
  dermal_cluster <- lapply(data_S_list, function(x){
    gene.ratio.cluster <- by(x@data["Col1a1",], INDICES = x@ident, FUN = function(xx) mean(xx > 0))
    return(names(gene.ratio.cluster)[gene.ratio.cluster > 0.85])
  })
  data_derm_S_list <- list()
  for(i in 1:n_sample){
    data_derm_S_list[[i]] <- SubsetData(
      data_S_list[[i]], cells.use = data_S_list[[i]]@cell.names[data_S_list[[i]]@ident %in% dermal_cluster[[i]]],
      subset.raw = T
    )
  }
  
  # Merge datasets
  E14_derm_S <- MergeSeurat(
    object1 = data_derm_S_list[[1]], object2 = data_derm_S_list[[2]]
  )
  E14_derm_S@var.genes <- union(data_derm_S_list[[1]]@var.genes, data_derm_S_list[[2]]@var.genes)
  E14_derm_S <- ScaleData(E14_derm_S)
  
  E14_derm_S@meta.data$wnt <- "WT"
  data_derm_S_list[[3]]@meta.data$wnt <- "LOF"
  
  data_derm_E14_S <- RunCCA(
    object = E14_derm_S, object2 = data_derm_S_list[[3]], num.cc = 15, group.by = "orig.ident"
  )
  data_derm_E14_S <- AlignSubspace(data_derm_E14_S, grouping.var = "orig.ident", dims.align = 1:15)
  
  # t-SNE on CCA.aligned
  data_derm_E14_S <- RunTSNE(
    data_derm_E14_S, reduction.use = "cca.aligned", dims.use = 1:15, reduction.name = "tsne", 
    check_duplicates = FALSE, do.fast = TRUE, seed.use=3, tsne.method="FIt-SNE", 
    fast_tsne_path="~/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
  )
  
  saveRDS(data_derm_E14_S, 'data/derm_S.rds')
  data_S <- data_derm_E14_S
  rm(data_derm_E14_S)
} else {
  data_S <- readRDS('data/derm_S.rds')
}



## Functions

# perform ALRA on Seurat object
alraSeurat2 <- function(obj, id = NULL, ...){
  if(is.null(id)){
    data_alra <- t(alra(t(as.matrix(obj@data)), ...)[[3]])
    colnames(data_alra) <- obj@cell.names
    data_alra <- Matrix(data_alra, sparse = T)
  } else {
    all.id <- unique(obj@meta.data[,id])
    data_alra_l <- list()
    # ALRA on each id
    for(ii in 1:length(all.id)){
      current.id <- all.id[ii]
      data_alra_l[[ii]] <- t(alra(t(as.matrix(obj@data[,obj@meta.data[,id] == current.id])), ...)[[3]])
      colnames(data_alra_l[[ii]]) <- obj@cell.names[obj@meta.data[,id] == current.id]
    }
    
    # combine data
    data_alra <- data_alra_l[[1]]
    for(ii in 2:length(all.id)){
      data_alra <- cbind(data_alra, data_alra_l[[ii]])
    }
    
    # reorder cells
    data_alra <- data_alra[,obj@cell.names]
    data_alra <- Matrix(data_alra, sparse = T)
  }
  
  obj@data <- data_alra
  return(obj)
}

# get Wnt pos ratio
getWnt <- function(obj){
  obj@meta.data$wntact <- "Other"
  obj@meta.data$wntact[obj@data["Axin2",] > 0 & obj@data["Lef1",] > 0] <- "Wnt Pos"
  
  if(is.null(obj@misc)){
    obj@misc <- list()
  }
  obj@misc[["Wnt"]] <- c(
    "WT" = sum(obj@meta.data$wntact == "Wnt Pos" & obj@meta.data$wnt == "WT") / sum(obj@meta.data$wnt == "WT"), 
    "Mutant" = sum(obj@meta.data$wntact == "Wnt Pos" & obj@meta.data$wnt == "LOF") / sum(obj@meta.data$wnt == "LOF")
  )
  return(obj)
}


# get DC cell ratio
getDC <- function(obj){
  if(is.null(obj@misc)){
    obj@misc <- list()
  }
  
  obj@misc[["DC"]] <- c(
    "WT" = sum(obj@data["Sox2",] > 0 & obj@meta.data$wnt == "WT") / sum(obj@meta.data$wnt == "WT"), 
    "Mutant" = sum(obj@data["Sox2",] > 0 & obj@meta.data$wnt == "LOF") / sum(obj@meta.data$wnt == "LOF")
  )
  return(obj)
}


# density plot for genes
cmpDensity <- function(datain, labelin, xlab = NULL, ylab = NULL, title = NULL){
  # plotidx <- which(datain != 0)
  # datain <- datain[plotidx]
  # labelin <- labelin[plotidx]
  myggdata <- data.frame(exp = datain, label = labelin)
  myggplot <- ggplot(myggdata) + geom_density(aes(exp, ..density.., fill = label), alpha = 0.5) + theme_cowplot() + 
    theme(legend.title = element_blank()) + xlab(xlab) + ylab(ylab) + ggtitle(title)
  return(myggplot)
}


# produce three density plots
get3density <- function(obj){
    t <- theme(legend.position = "none", axis.text = element_text(size=10), plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"))
  mygg1 <- cmpDensity(datain = obj@data["Axin2",], labelin = obj@meta.data$wnt,title= 'Axin2') + t
  mygg2 <- cmpDensity(datain = obj@data["Lef1",], labelin = obj@meta.data$wnt, title='Lef1') + t
  mygg3 <- cmpDensity(datain = obj@data["Sox2",], labelin = obj@meta.data$wnt, title='Sox2') + t + theme(legend.position = c(0.6,0.8), legend.text=element_text(size=10), legend.spacing.x = unit(0.1,'cm'))
  #mygg <- grid.arrange(grobs = list(mygg1,mygg2,mygg3), nrow = 1)
  mygg <- list(mygg1,mygg2,mygg3)
  return(mygg)
}



## Original

# data_S <- readRDS("./data_derm_E14_S.rds")
striptitle <- function(x)  {
    plotLabels <- x$labels
    plotLabels$title <- NULL
    x$labels <- plotLabels 
    x
}

TSNEPlot(data_S, group.by = "wnt", pt.size = 0.5)

gg1 <- lapply(FeaturePlot(data_S, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                   pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                     x + theme_cowplot() + 
                       theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                             axis.ticks = element_blank())
                   })
gg2 <- get3density(data_S)
g <-  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2)
plot(g)
ggsave(
 g, 
  filename = "./figs/original_combine.pdf", width = 6, height = 4
)
plot(striptitle(gg1$Axin2))
lapply(gg1,striptitle)

plot(gg1$Axin2)

# plot legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
ggsave(g_legend(gg2[[1]] + theme(legend.position = "top")), filename = "./figs/density_legend.pdf", width = 1.2, height = 0.5)

data_S <- getWnt(data_S)
TSNEPlot(data_S, group.by = "wntact", colors.use = c("gray","darkgreen"), pt.size = 0.5)
data_S <- getDC(data_S)



## Imputation method

# ALRA
data_S_alra <- alraSeurat2(data_S, id = "wnt")

gg1 <- lapply(FeaturePlot(data_S_alra, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                   pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                     x + theme_cowplot() + 
                       theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                             axis.ticks = element_blank())
                   })
gg2 <- get3density(data_S_alra)
ggsave(
  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2), 
  filename = "./figs/alra_combine.pdf", width = 6, height = 4
)

data_S_alra <- getWnt(data_S_alra)
data_S_alra <- getDC(data_S_alra)



# MAGIC
data_norm_WT_magic <- t(magic(t(as.matrix(data_S@data[,data_S@meta.data$wnt == "WT"])))$result)
data_norm_LOF_magic <- t(magic(t(as.matrix(data_S@data[,data_S@meta.data$wnt == "LOF"])))$result)

data_norm_magic <- cbind(data_norm_WT_magic, data_norm_LOF_magic)[,data_S@cell.names]

# data_norm_magic <- t(magic(t(as.matrix(data_S@data)))$result)

data_S_magic <- data_S
data_S_magic@data <- Matrix(data_norm_magic)

gg1 <- lapply(FeaturePlot(data_S_magic, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                   pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                     x + theme_cowplot() + 
                       theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                             axis.ticks = element_blank())
                   })
gg2 <- get3density(data_S_magic)
ggsave(
  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2), 
  filename = "./figs/magic_combine.pdf", width = 6, height = 4
)

data_S_magic <- getWnt(data_S_magic)
data_S_magic <- getDC(data_S_magic)



# DCA
fn =  "data/data_filtered_WT.csv"
if ( !file.exists(fn)) {
  write.table( as.matrix(data_S@raw.data[rownames(data_S@data),data_S@meta.data$wnt == "WT"]), fn, 
               sep = ",", quote = F, row.names = T, col.names = T
  )
}else {
  print("Using existing filtered CSV instead of writing new one.")
}

fno =  "data/dca_out_WT/mean.tsv"
if ( !file.exists(fno)) {
  print(sprintf("Generating %s", fno))
  system('dca data/data_filtered_WT.csv data/dca_out_WT')
}
data_norm_WT_dca <- as.matrix(read.table(fno, sep = "\t", header = T, row.names = 1, stringsAsFactors = F))
data_norm_WT_dca <- log(data_norm_WT_dca + 1)
colnames(data_norm_WT_dca) <- gsub(".","-",colnames(data_norm_WT_dca), fixed = T)

fn =  "data/data_filtered_LOF.csv"
if ( !file.exists(fn)) {
  write.table( as.matrix(data_S@raw.data[rownames(data_S@data),data_S@meta.data$wnt == "LOF"]), fn, 
               sep = ",", quote = F, row.names = T, col.names = T
  )
}else {
  print("Using existing filtered CSV instead of writing new one.")
}

fno =  "data/dca_out_LOF/mean.tsv"
if ( !file.exists(fno)) {
  print(sprintf("Generating %s", fno))
  system('dca data/data_filtered_LOF.csv data/dca_out_LOF')
}
data_norm_LOF_dca <- as.matrix(read.table(fno, sep = "\t", header = T, row.names = 1, stringsAsFactors = F))
data_norm_LOF_dca <- log(data_norm_LOF_dca + 1)
colnames(data_norm_LOF_dca) <- gsub(".","-",colnames(data_norm_LOF_dca), fixed = T)

# function to combine data
combineData <- function(data1, data2, name1 = NULL, name2 = NULL, delim = "_"){
  exp_pre <- data1
  exp_post <- data2
  if(!is.null(name1)){
    colnames(exp_pre) <- paste(colnames(exp_pre), name1, sep = delim)
  }
  if(!is.null(name2)){
    colnames(exp_post) <- paste(colnames(exp_post), name2, sep = delim)
  }
  
  # find all genes
  allgene <- unique(c(rownames(exp_pre), rownames(exp_post)))
  missG_pre <- setdiff(allgene, rownames(exp_pre))
  missG_post <- setdiff(allgene, rownames(exp_post))
  
  # create matrix with 0 for missing genes
  miss_pre <- matrix(0, ncol = dim(exp_pre)[2], nrow = length(missG_pre))
  rownames(miss_pre) <- missG_pre
  colnames(miss_pre) <- colnames(exp_pre)
  miss_post <- matrix(0, ncol = dim(exp_post)[2], nrow = length(missG_post))
  rownames(miss_post) <- missG_post
  colnames(miss_post) <- colnames(exp_post)
  
  # bind data
  new_pre <- rbind(exp_pre, miss_pre)
  new_post <- rbind(exp_post, miss_post)
  new_pre <- new_pre[order(rownames(new_pre)),]
  new_post <- new_post[order(rownames(new_post)),]
  print(mean(rownames(new_pre) == rownames(new_post)))
  data_new <- cbind(new_pre, new_post)
  return(data_new)
}


data_norm_dca <- combineData(data1 = data_norm_WT_dca, data_norm_LOF_dca)[,data_S@cell.names]
data_S_dca <- data_S
data_S_dca@raw.data <- data_norm_dca
data_S_dca@data <- data_norm_dca
data_S_dca <- NormalizeData(data_S_dca)

gg1 <- lapply(FeaturePlot(data_S_dca, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                   pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                     x + theme_cowplot() + 
                       theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                             axis.ticks = element_blank())
                   })
gg2 <- get3density(data_S_dca)
ggsave(
  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2), 
  filename = "./figs/dca_combine.pdf", width = 6, height = 4
)

data_S_dca <- getWnt(data_S_dca)
data_S_dca <- getDC(data_S_dca)



## SAVER
set.seed(3)
saver_WT <- saver(
  data_S@raw.data[rownames(data_S@data),data_S@meta.data$wnt == "WT"],
  pred.genes = match(c("Axin2","Lef1","Sox2"), rownames(data_S@data)), pred.genes.only = T
)
saver_norm_WT <- sample.saver(saver_WT, seed = 3)

saver_LOF <- saver(
  data_S@raw.data[rownames(data_S@data),data_S@meta.data$wnt == "LOF"],
  pred.genes = match(c("Axin2","Lef1","Sox2"), rownames(data_S@data)), pred.genes.only = T
)
saver_norm_LOF <- sample.saver(saver_LOF, seed = 3)

saver_norm <- log(cbind(saver_norm_WT, saver_norm_LOF)[,data_S@cell.names] + 1)
data_S_saver <- data_S
data_S_saver@data <- saver_norm

gg1 <- lapply(FeaturePlot(data_S_saver, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                          pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                            x + theme_cowplot() + 
                              theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                                    axis.ticks = element_blank())
                          })
gg2 <- get3density(data_S_saver)
ggsave(
  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2), 
  filename = "./figs/saver_combine.pdf", width = 6, height = 4
)

data_S_saver <- getWnt(data_S_saver)
data_S_saver <- getDC(data_S_saver)



## scImpute

scimpute( count_path = "data/data_filtered_WT.csv",
          out_dir = "data/scimpute_out_WT/", 
          Kcluster = 3)
scimpute_WT <- as.matrix(read.table( "data/scimpute_out_WT/scimpute_count.csv", 
                                       sep = ",", header = T, row.names = 1,
                                       stringsAsFactors = F))
colnames(scimpute_WT) <- gsub("\\.","-", colnames(scimpute_WT))

scimpute( count_path = "data/data_filtered_LOF.csv",
          out_dir = "data/scimpute_out_LOF/", 
          Kcluster = 3)
scimpute_LOF <- as.matrix(read.table( "data/scimpute_out_LOF/scimpute_count.csv", 
                                     sep = ",", header = T, row.names = 1,
                                     stringsAsFactors = F))
colnames(scimpute_LOF) <- gsub("\\.","-", colnames(scimpute_LOF))

all(rownames(scimpute_WT) == rownames(scimpute_LOF))
scimpute_norm <- cbind(scimpute_WT,scimpute_LOF)[,data_S@cell.names]

data_S_scimpute <- data_S
data_S_scimpute@raw.data <- scimpute_norm
data_S_scimpute@data <- scimpute_norm
data_S_scimpute <- NormalizeData(data_S_scimpute)

gg1 <- lapply(FeaturePlot(data_S_scimpute, features.plot = c("Axin2","Lef1","Sox2"), cols.use = c("gray","red"), 
                          pt.size = 0.5, nCol = 3, no.axes = T, no.legend = T, do.return = T), function(x) {
                            x + theme_cowplot() + 
                              theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(),
                                    axis.ticks = element_blank())
                          })
gg2 <- get3density(data_S_scimpute)
ggsave(
  grid.arrange(grobs = c( gg2, lapply(gg1, striptitle)), nrow = 2), 
  filename = "./figs/scimpute_combine.pdf", width = 6, height = 4
)

data_S_scimpute <- getWnt(data_S_scimpute)
data_S_scimpute <- getDC(data_S_scimpute)



## Print results

print(data_S_alra@misc)
print(data_S_magic@misc)
print(data_S_dca@misc)
print(data_S_saver@misc)
print(data_S_scimpute@misc)
