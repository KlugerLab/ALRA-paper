### Compare imputation methods on Manolis ds005 data

library(ggplot2)
library(cowplot)
library(gridExtra)
library(Seurat) #V2.3.4

source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")


##=============================================##
## Data prep


## Data download

# download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185638
# and save to data/ folder


## Load data

data_exp <- read.table("data/GSM5621018_Run5C_gene_exon_tagged_dge.txt.gz",
                       sep = "\t", header = T, as.is = T, row.names = 1)


# Seurat
data_S <- CreateSeuratObject(raw.data = data_exp, min.genes = 0)

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
data_original_S <- data_S




##=============================================##
## Imputation

## alra

fn ='data/ds005_alra.RDS';
if ( !file.exists(fn)) {
  set.seed(3)
  chose.k <- choose_k(as.matrix(t(data_original_S@data)))
  saveRDS(chose.k,'data/ds005_alra_choose_k.RDS')
  data_alra_S <- alraSeurat2(data_original_S,k=chose.k$k)
  saveRDS(data_alra_S,fn)
}else{
  print(sprintf("Loading %s\n", fn));
  chose.k <- readRDS('data/ds005_alra_choose_k.RDS')
  data_alra_S <- readRDS(fn)
}



## magic
library(Rmagic)

fn ='data/ds005_magic.RDS';
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
fn <- "data/dca_ds005/mean.tsv"
if (!file.exists(fn)) {
  write.table(
    data_original_S@raw.data[rownames(data_original_S@data),data_original_S@cell.names],
    file = "data/ds005_count.txt", 
    sep = "\t", quote = F, row.names = T, col.names = T
  )
  system.time(system(sprintf('dca data/ds005_count.txt data/dca_ds005')))
}

data_dca <- read.table( fn, sep="\t", header = TRUE, row.names = 1, stringsAsFactors = F)
data_dca <- Matrix(as.matrix(data_dca[rownames(data_original_S@data),data_original_S@cell.names]))
rownames(data_dca) <- rownames(data_original_S@data)
data_dca[is.na(data_dca)] <- 0
data_dca_S <- CreateSeuratObject(raw.data = data_dca)
data_dca_S <- NormalizeData(data_dca_S)



## SAVER

library(SAVER)
fn ='data/ds005_saver.rds';
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
fn ='data/ds005_scimpute.rds';
if ( !file.exists(fn)) {
  cat(sprintf("Generating %s\n", fn))
  write.table( data_original_S@raw.data[rownames(data_original_S@data),data_original_S@cell.names], 
               "data/ds005_count.csv", sep = ",", quote = F, row.names = T, col.names = T)
  scimpute( count_path = "data/ds005_count.csv",
            out_dir = "data/scimpute_ds005/", 
            Kcluster = 10)
  scimpute_norm <- as.matrix(read.table( "data/scimpute_ds005/scimpute_count.csv", 
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
    geom_density( fill="#3182bd",data = data.frame(x = x[x>0]), aes(x)) + 
    xlab(xlab)
  return(myggplot)
}
# plot Pdgfra
g <- grid.arrange(grobs = list(
  plotDensity(data_original_S@data["Pdgfra",], "Original"),
  plotDensity(data_alra_S@data["Pdgfra",], "ALRA"),
  plotDensity(data_magic_S@data["Pdgfra",], "MAGIC"),
  plotDensity(data_dca_S@data["Pdgfra",], "DCA"),
  plotDensity(data_saver_S@data["Pdgfra",], "SAVER"),
  plotDensity(data_scimpute_S@data["Pdgfra",], "scImpute")
), ncol = 3)
ggsave(g, width=12, height=6,filename = "figs/methods_pdgfra_ds005.pdf")


## Flow result
flow_pdgfra <- read.csv('data/PdgfraEGFP_histogram_data.csv', stringsAsFactors=F)
flow_pdgfra2 <- rep(flow_pdgfra$Compensated.Al488.channel, times = flow_pdgfra$Number.of.cells)

myggplot <- ggplot() + theme_cowplot() + 
  geom_density( fill="#3182bd",data = data.frame(x = flow_pdgfra2), aes(x)) + scale_x_log10() + xlab ('EGFP fluorescence')
ggsave(myggplot, filename = "figs/flow_cytometry_pdgfra.pdf", width = 6, height = 4)


