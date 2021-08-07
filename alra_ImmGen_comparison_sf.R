library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(cowplot)
library(rsvd)
source("../ALRA/alra.R")
source("alra_ImmGen_funs.R")
##=======================================================##
## Generate simulation data


## Get p for each cell type from bulk RNA-seq

# load ImmGen bulk RNA-seq data
load("data.bu/ImmGen_bulk.RData")
colnames(bulkCount) <- bulkType
nCell <- 5000
nGene <- nrow(bulkCount)

# load cell size data
load("data.bu/PBMC_lib_size.RData")
ss <- which(!bulkType %in%   c('B.Cell', 'NK.Cell', 'Neutrophil', 'Dendritic.Cell') )
#ss <-  which(!bulkType %in% c('fofo'))
bulkType <- bulkType[ss]
bulkCount  <- bulkCount[,ss]
bulkNorm  <- bulkNorm[,ss]
set.seed(3)
# subsample from PBMC library size
simCellSize <- sample(PBMC_lib_size, size = nCell, replace = F)

# get prob
bulkProb <- sweep(bulkCount, MARGIN = 2, STATS = colSums(bulkCount), FUN = "/")
cellType <- sample(bulkType, size = nCell, replace = T)

# simulate cell type
true.mat <- bulkCount[,cellType]
true.p <- sweep(true.mat,2, colSums(true.mat), '/')

true.mat[1:5,1:5]

set.seed(3)
frac <- 1
simCellSize_ <-frac*simCellSize
simCount <- matrix(0, nrow = nrow(true.p), ncol = ncol(true.p))
for(i in 1:ncol(true.p)){
  simCount[,i] <- rmultinom(n = 1, size = simCellSize_[i], prob = true.p[,i])
}
keep.cells <- which ( colSums(simCount) > 0)
cellType.sub <- cellType[keep.cells]
keep.genes <- which ( rowSums(simCount) > 0)
simCount <- simCount[keep.genes,keep.cells]
true.p.sub <- true.p[keep.genes,keep.cells]
dim(simCount)
sum(rowSums(simCount >0) == 0)

sfs <- c(1000,5000,10000,20000)
true.nonzeros <- sfs
true.zeros <- sfs
i <- 1
for(sf in sfs){
  simNorm <- log(normbycol(simCount, newsum = sf) + 1)
  cat(sprintf("Max = %s\n", max(simNorm)))
  ## ALRA
  set.seed(3)
  simAlraK    <- choose_k(t(as.matrix(simNorm)), q = 10)
  cat(sprintf("Choose k = %s\n", simAlraK$k))
  simAlra     <- alra(t(as.matrix(simNorm)), k = simAlraK$k)
  simNormAlra <- t(simAlra$A_norm_rank_k_cor_sc)
  
  true.nonzeros[i]  <- sum(simNormAlra >0  &  true.p.sub >0)/sum(true.p.sub>0)
  true.zeros[i]  <- sum(simNormAlra <=0  &  true.p.sub ==0)/sum(true.p.sub==0)
  i <- i+1
}

print(format(true.zeros, digits=3))
print(format(true.nonzeros, digits=3))
