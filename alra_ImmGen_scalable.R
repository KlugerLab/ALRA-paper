### alra_ImmGen_scalable.R
### @Jun Zhao
### script to generate scalable results of ALRA with different cell number

library(ggplot2)
library(reshape2)
library(rsvd)

source("./alra.R")
source("./alra_ImmGen_funs.R")


##=======================================================##
## Generate simulation data

set.seed(3)

## Get p for each cell type from bulk RNA-seq

# load ImmGen bulk RNA-seq data
load("./ImmGen_bulk.RData")

# get prob
bulkProb <- sweep(bulkCount, MARGIN = 2, STATS = colSums(bulkCount), FUN = "/")

# simulation data init
nCell <- 60000
nGene <- nrow(bulkCount)



## Get cell transcript count from PBMC data

# load cell size data
load("./PBMC_lib_size.RData")

# subsample from PBMC library size
simCellSize <- sample(PBMC_lib_size, size = nCell, replace = F)



## Simulation

# simulation data matrix
simCount <- matrix(0, nrow = nGene, ncol = nCell)

# simulate cell type
cellType <- sample(bulkType, size = nCell, replace = T)
table(cellType)


# get count from multinomial distribution
for(i in 1:nCell){
  simCount[,i] <- rmultinom(n = 1, size = simCellSize[i], prob = bulkProb[,which(bulkType == cellType[i])])
}

sum(colSums(simCount) >= 1000) # 50832

rownames(simCount) <- rownames(bulkCount)
colnames(simCount) <- paste("Cell", seq(1:nCell), sep = "_")

simCount <- simCount[rowSums(simCount) > 0,]
mean(simCount == 0) # 0.9427712





##=======================================================##
## 

## Settings

nRealCell <- c(1000,5000,10000,20000,30000,40000,50000)
nReal <- length(nRealCell)

zero.eval.list.ALRA <- list()
zero.eval.list.observed <- list()



## Go through all cell number

for(i in 1:nReal){
  set.seed(3)
  
  realCellIdx <- sample(which(colSums(simCount) >= 1000), size = nRealCell[i], replace = F)
  print(mean(simCount[,realCellIdx] == 0))
  simNorm <- log(normbycol(simCount[rowSums(simCount[,realCellIdx]) != 0,realCellIdx], 10000) + 1)
  
  simAlraK <- choose_k(t(as.matrix(simNorm)))
  simAlra <- alra(t(as.matrix(simNorm)), k = simAlraK$k)
  simNormAlra <- t(simAlra$A_norm_rank_k_cor_sc)
  
  zero.eval.list.ALRA[[as.character(nRealCell[i])]] <- zeroEval(true.data = bulkCount[rownames(simNorm),], 
                                                                true.label = bulkType, 
                                                                impute.data = simNormAlra, 
                                                                impute.label = cellType[realCellIdx])
  
  zero.eval.list.observed[[as.character(nRealCell[i])]] <- zeroEval(true.data = bulkCount[rownames(simNorm),], 
                                                                    true.label = bulkType, 
                                                                    impute.data = simNorm, 
                                                                    impute.label = cellType[realCellIdx])
  
  # subtract observed non zero ratio
  zero.eval.list.ALRA[[as.character(nRealCell[i])]][,2] <- zero.eval.list.ALRA[[as.character(nRealCell[i])]][,2] - 
    zero.eval.list.observed[[as.character(nRealCell[i])]][,2]
}



## Plot

# zero.ratio
zero.ratio.mat.ALRA <- matrix(unlist(lapply(zero.eval.list.ALRA, FUN = function(x){return(x[,1])})), 
                              nrow = length(zero.eval.list.ALRA), byrow = T)
colnames(zero.ratio.mat.ALRA) <- rownames(zero.eval.list.ALRA$`1000`)
rownames(zero.ratio.mat.ALRA) <- names(zero.eval.list.ALRA)

gg1 <- ggplot(data = melt(zero.ratio.mat.ALRA), 
              aes(x = Var1, y = 100*(1-value), group = Var2, col = factor(as.character(Var2)))) + 
  theme_cowplot() + 
  geom_line() + scale_x_continuous(labels = c("1k","5k","10k","20k","30k","40k","50k"), breaks = nRealCell) + 
  geom_point() + xlab("N") + ylab("Biological zeros completed (%)") + 
  theme(legend.position = "none")


# non.zero.ratio
non.zero.ratio.mat.ALRA <- matrix(unlist(lapply(zero.eval.list.ALRA, FUN = function(x){return(x[,2])})), 
                                  nrow = length(zero.eval.list.ALRA), byrow = T)
colnames(non.zero.ratio.mat.ALRA) <- rownames(zero.eval.list.ALRA$`1000`)
rownames(non.zero.ratio.mat.ALRA) <- names(zero.eval.list.ALRA)

gg2 <- ggplot(data = melt(non.zero.ratio.mat.ALRA), 
              aes(x = Var1, y = 100*value, group = factor(Var2), 
                  col = factor(as.character(Var2)))) + 
  theme_cowplot() + 
  geom_line() + scale_x_continuous(labels = c("1k","5k","10k","20k","30k","40k","50k"), breaks = nRealCell) + 
  geom_point() + xlab("N") + ylab("Technical zeros completed (%)") + 
  theme(legend.position = "none")

ggsave(gridExtra::grid.arrange(gg1,gg2,nrow = 1), filename = "./ImmGen_scalable_zero.pdf", 
       width = 8, height = 4)


# make legend
gg3 <- ggplot(data = melt(non.zero.ratio.mat.ALRA), 
              aes(x = Var1, y = 100*value, 
                  group = factor(Var2), 
                  col = factor(gsub("."," ",as.character(Var2),fixed = T)))) + theme_cowplot() + 
  geom_line() + scale_x_discrete(labels = c("0"="2,500","1"="5,000","2"="10,000","3"="20,000","4"="40,000")) + 
  geom_point() + xlab("N") + ylab("Technical zeros completed (%)") + 
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(grid::grid.draw(g_legend(gg3)), filename = "./ImmGen_zero_legend.pdf", width = 8, height = 2)



