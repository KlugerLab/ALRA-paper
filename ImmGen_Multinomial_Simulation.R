library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(reshape2)
library(cowplot)
library(rsvd)
source("../ALRA/alra.R")
source("alra_ImmGen_funs.R")
##=======================================================##
## Generate simulation data

generateSimulatedData  <- function (  nCell = 5000, sf = 10000 ) {
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
    simNorm <- log(normbycol(simCount, newsum = sf) + 1)
    list( simNorm = simNorm, true.p = true.p.sub)
}

simulatedData  <- generateSimulatedData()
A <-  simulatedData$simNorm
P <-  simulatedData$true.p
k    <- choose_k(t(as.matrix(A)), q = 10)
A.alra     <- alra(t(as.matrix(A)), k = k$k)[[3]]
dim(A.alra)
imputed.values  <- data.frame(P=c(P), A=c(A), A.alra=c(t(A.alra)), A.norm = log(10000*c(A) + 1), P.norm = log(10000*c(P) + 1) ) 
imputed.values  <- imputed.values  %>% mutate( to.impute = P >0 & A == 0, to.preserve = P==0  )
imputed.values.filtered <- imputed.values %>% filter (P >0 & A == 0) 
g <- ggplot(imputed.values.filtered, aes(x=P, y=1.0*(A.alra>0))) + geom_smooth() + scale_x_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3)) + cowplot::theme_cowplot() + labs (x='Probability of expression', y='Probability of imputation') + ylim(0,1)
ggsave(g, file = 'figs/prob.recovered.pdf')




#imputed.values$P.quantile  <- cut(imputed.values$P, breaks = quantile( imputed.values$P, probs=seq(0,1,by = 0.1))) imputed.values %>% group_by(P.quantile) %>% summarise(mean(A.alra), mean(A.alra>0), mean(P.norm) ) #with(imputed.values,  mean( A.alra[to.preserve] >0 )) with(imputed.values, cor( P.norm[P.norm >0], A.alra[P.norm >0])) dim(A) with(imputed.values, sum(P.norm[P.norm>0] >0 & P.norm[ P[1:50,1:50] A.alra[1:50,1:50]

