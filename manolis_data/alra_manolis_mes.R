### Compare imputation methods on colon mesenchyme data from Roulis et al. (2020) Nature

library(ggplot2)
library(cowplot)
library(gridExtra)
library(Seurat) #V2.3.4

source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")



##=============================================##
## Functions

# Combine data from multiple files
combineData <- function(data1, data2, name1 = NULL, name2 = NULL, delim = "_"){
  # read in data
  # exp_pre <- read.table(data1, header = T,sep = "\t",stringsAsFactors = F, row.names = 1)
  # exp_post <- read.table(data2, header = T,sep = "\t",stringsAsFactors = F, row.names = 1)
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





##=============================================##
## Data prep

## Download data

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE142423&format=file", 
  destfile = "./data/GSE142423_RAW.tar"
)
system(
  "tar -xvf ./data/GSE142423_RAW.tar -C ./data/"
)


## Prepare data

# samples
sample_names <- c("C1A","C1B","C2A","C2B","C2C")
n_sample <- length(sample_names)

# read in data
data_list <- vector("list", n_sample)
names(data_list) <- sample_names

for(i in 1:n_sample){
  data_list[[i]] <- read.table(
    grep(sample_names[i], list.files("./data/", full.names = T), value = T), 
    sep = "\t", row.names = 1, header = T, stringsAsFactors = F
  )
}



## Combine data directly

data_exp <- combineData(
  data1 = data_list[[1]], data2 = data_list[[2]], name1 = sample_names[1], name2 = sample_names[2]
)

for(i in 3:n_sample){
  data_exp <- combineData(
    data1 = data_exp, data2 = data_list[[i]], name2 = sample_names[i]
  )
}


# Seurat
data_S <- CreateSeuratObject(raw.data = data_exp, min.genes = 0, names.delim = "_", names.field = 2)


# transcript sum filtering
sum(data_S@meta.data$nUMI >= 1000)
data_S <- FilterCells(data_S, subset.names = "nUMI", low.thresholds = 999)

table(data_S@meta.data$orig.ident)
head(data_S@meta.data)
min(data_S@meta.data$nUMI)


# filter cells on mito.ratio
mito.genes <- grep(pattern = "mt-", x = rownames(x = data_S@data), value = TRUE, ignore.case = T)
data_S@meta.data$mito.ratio <- colSums(data_S@raw.data[mito.genes,data_S@cell.names])/
  colSums(data_S@raw.data[,data_S@cell.names])
max(data_S@meta.data$mito.ratio)
VlnPlot(data_S, "mito.ratio", group.by = "orig.ident", point.size.use = 0.1)

data_S <- FilterCells(data_S, subset.names = "mito.ratio", high.thresholds = 0.1)


# analysis
data_S <- NormalizeData(data_S)
data_S <- ScaleData(data_S, do.scale = F, do.center = T)

data_S <- RunPCA(data_S, pcs.compute = 36, pc.genes = rownames(data_S@data), do.print = F)
plot(data_S@dr$pca@sdev)

data_S <- RunTSNE(
  data_S, dims.use = 1:36, seed.use = 1
)
TSNEPlot(data_S, group.by = "orig.ident", pt.size = 0.5)

data_original_S <- data_S




##=============================================##
## Imputation

## alra

fn ='data/mes_alra.RDS';
if ( !file.exists(fn)) {
  set.seed(3)
  chose.k <- choose_k(as.matrix(t(data_original_S@data)))
  saveRDS(chose.k,'data/mes_alra_choose_k.RDS')
  data_alra_S <- alraSeurat2(data_original_S,k=chose.k$k)
  saveRDS(data_alra_S,fn)
}else{
  print(sprintf("Loading %s\n", fn));
  chose.k <- readRDS('data/mes_alra_choose_k.RDS')
  data_alra_S <- readRDS(fn)
}



## magic
library(Rmagic)

fn ='data/mes_magic.RDS';
if ( !file.exists(fn)) {
  set.seed(3)
  data_magic <- magic(t(data_original_S@data))
  data_magic_S <- data_original_S
  data_magic_S@data <- Matrix(as.matrix(t(data_magic$result)))
  saveRDS(data_magic_S,fn)
}else{
  print(sprintf("Loading %s\n", fn));
  data_magic_S <- readRDS(fn)
}



## DCA

# read in dca result
fn <- "data/dca_mes/mean.tsv"
if (!file.exists(fn)) {
  write.table(
    data_original_S@raw.data[rownames(data_original_S@data),data_original_S@cell.names],
    file = "data/mes_count.txt", 
    sep = "\t", quote = F, row.names = T, col.names = T
  )
  system.time(system(sprintf('dca data/mes_count.txt data/dca_mes')))
}

data_dca <- read.table( fn, sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F)
data_dca <- Matrix(as.matrix(data_dca[rownames(data_original_S@data),data_original_S@cell.names]))
rownames(data_dca) <- rownames(data_original_S@data)
data_dca[is.na(data_dca)] <- 0
data_dca_S <- CreateSeuratObject(raw.data = data_dca)
data_dca_S <- NormalizeData(data_dca_S)



## SAVER

library(SAVER)
fn ='data/mes_saver.rds';
if ( !file.exists(fn)) {
  cat(sprintf("Generating %s\n", fn))
  saver_result <- saver(
    data_original_S@raw.data[rownames(data_original_S@data),data_original_S@cell.names]
  )
  saveRDS(saver_result,fn)
}else{
  print(sprintf("Loading %s\n", fn));
  saver_result <- readRDS(fn)
}
saver_exp <- sample.saver(x = saver_result, seed = 3)
saver_exp[is.na(saver_exp)] <- 0
data_saver_S <- CreateSeuratObject(
  raw.data = saver_exp, names.delim = "_", names.field = 2
)
data_saver_S <- NormalizeData(data_saver_S)
rm(saver_result)



## scImpute

library(scImpute)
fn ='data/mes_scimpute.rds';
if ( !file.exists(fn)) {
  cat(sprintf("Generating %s\n", fn))
  if(!file.exists("data/mes_count.csv")){
    write.table(
      data_original_S@raw.data[rownames(data_original_S@data),data_original_S@cell.names],
      file = "data/mes_count.csv", 
      sep = ",", quote = F, row.names = T, col.names = T
    )
  }
  scimpute( count_path = "data/mes_count.csv",
            out_dir = "data/scimpute_mes/", 
            Kcluster = 16)
  scimpute_norm <- as.matrix(read.table( "data/scimpute_mes/scimpute_count.csv", 
                                         sep = ",", header = T, row.names = 1,
                                         stringsAsFactors = F))
  rownames(scimpute_norm) <- gsub("_","-", rownames(scimpute_norm))
  data_scimpute_S <- CreateSeuratObject(
    raw.data = scimpute_norm, names.delim = "_", names.field = 2
  )
  data_scimpute_S <- NormalizeData(data_scimpute_S)
  saveRDS(data_scimpute_S,fn)
}else{
  print(sprintf("Loading %s\n", fn));
  data_scimpute_S <- readRDS(fn)
}





##=============================================##
## Pdgfra

# function
plotDensity <- function(x, xlab){
  myggplot <- ggplot() + theme_cowplot() + 
    geom_density(fill="#3182bd",data = data.frame(x = x[x>0]), aes(x)) + 
    xlab(xlab)
  return(myggplot)
}
# plot Pdgfra
g<- grid.arrange(grobs = list(
  plotDensity(data_original_S@data["Pdgfra",], "Original"),
  plotDensity(data_alra_S@data["Pdgfra",], "ALRA"),
  plotDensity(data_magic_S@data["Pdgfra",], "MAGIC"),
  plotDensity(data_dca_S@data["Pdgfra",], "DCA"),
  plotDensity(data_saver_S@data["Pdgfra",], "SAVER"),
  plotDensity(data_scimpute_S@data["Pdgfra",], "scImpute")
), ncol = 3)
ggsave(g, width=12, height=6,filename = "figs/methods_pdgfra_mes.pdf")


