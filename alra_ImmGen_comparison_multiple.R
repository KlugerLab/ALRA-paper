library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(cowplot)
library(rsvd)
library(Rmagic)
library(SAVER)
library(scImpute)
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



fracs  <- c(0.25, 0.5,0.75, 1)
read.depths <- fracs
methodslist <- c('ALRA','MAGIC', 'SAVER', 'scImpute')
true.zeros <- setNames( data.frame(matrix(-1, ncol=2, nrow = length(fracs))), methodslist)
true.nonzeros <- setNames(data.frame(matrix(-1, ncol=2, nrow = length(fracs))), methodslist)
set.seed(3)
for (fraci in 1:length(fracs)) {
    frac <- fracs[fraci]
    simCellSize_ <-frac*simCellSize
    read.depths[fraci] <- mean(simCellSize_)
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
    simNorm <- log(normbycol(simCount, 10000) + 1)
    ################################
    # Imputation  
    ################################
    ## ALRA
    simAlraK    <- choose_k(t(as.matrix(simNorm)))
    simAlra     <- alra(t(as.matrix(simNorm)), k = simAlraK$k)
    simNormAlra <- t(simAlra$A_norm_rank_k_cor_sc)
    #MAGIC
    reticulate::use_condaenv('magic')
    simMagic     <- magic(data = t(as.matrix(simNorm)))
    simNormMagic <- t(simMagic$result)
    #SAVER
    simSaver <- saver(x = simCount,ncores = 10 )
    simNormSaver <- sample.saver(simSaver)
    #scImpute
    write.table(
      simCount, "data/simulation_5000.csv", 
      sep = ",", quote = F, row.names = T, col.names = T
    )
    system('Rscript runscimpute.R')
    scimpute_norm <- as.matrix(read.table(
      "data/scimpute_simulation/scimpute_count.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F
    ))
    simNormScimpute <- log(normbycol(scimpute_norm, 10000) + 1)
    ################################
    # Evaluate 
    ################################
    true.nonzeros[fraci,1]  <- sum(simNormAlra >0  &  true.p.sub >0)/sum(true.p.sub>0)
    true.zeros[fraci,1]  <- sum(simNormAlra <=0  &  true.p.sub ==0)/sum(true.p.sub==0)
    true.nonzeros[fraci,2]  <- sum(simNormMagic >0  &  true.p.sub >0)/sum(true.p.sub>0)
    true.zeros[fraci,2]  <- sum(simNormMagic <=0  &  true.p.sub ==0)/sum(true.p.sub==0)
    dim(simNormSaver)
    dim(true.p.sub)
    true.nonzeros[fraci,3]  <- sum(simNormSaver >0  &  true.p.sub >0)/sum(true.p.sub>0)
    true.zeros[fraci,3]  <- sum(simNormSaver <=0  &  true.p.sub ==0)/sum(true.p.sub==0)
    true.nonzeros[fraci,4]  <- sum(simNormScimpute >0  &  true.p.sub >0)/sum(true.p.sub>0)
    true.zeros[fraci,4]  <- sum(simNormScimpute <=0  &  true.p.sub ==0)/sum(true.p.sub==0)
    print( true.zeros)
    print(true.nonzeros)
}
saveRDS(true.zeros,file = 'true.zeros.RDS')
saveRDS(true.nonzeros,file = 'true.nonzeros.RDS')

true.zeros.long <- cbind( fracs=read.depths, true.zeros) %>% pivot_longer( ! fracs, names_to = 'method', values_to = 'value')
g1 <- ggplot(true.zeros.long, aes(x=fracs, y=value, color=method)) + geom_point()  + geom_line() + ylab('Biological zeros preserved') + cowplot::theme_cowplot()+ theme(text = element_text(size = 11), axis.text = element_text(size=11))+ ylim(0,1)
true.nonzeros.long <- cbind( fracs=read.depths, true.nonzeros) %>% pivot_longer( ! fracs, names_to = 'method', values_to = 'value')
g2 <- ggplot(true.nonzeros.long, aes(x=fracs, y=value, color=method)) + geom_point()  + geom_line()+ ylab('Technical zeros imputed')  + cowplot::theme_cowplot() + theme(text = element_text(size = 11), axis.text = element_text(size=11)) + ylim(0,1)
g <- g1+g2 + plot_layout(guides = "collect")  & theme(legend.position = "bottom")& xlab('Average read depth')  
ggsave(g, file="figs/simulation_accuracy.pdf", height = 3)

#ggsave(g, file="figs/test.pdf")
#
#
#zero.eval.MAGIC
#mean(simNormMagic >0)
#
#toplot1 <- data.frame( y = true.zeros, x = 
#
#################################
## Look at it 
#################################
#row_gene_means <- rowMeans(simNorm)
#simNorm_c <- sweep(simNorm, 1, row_gene_means, '-')
#system.time(svdout <- rsvd(simNorm_c, k=100))
#obsPCA <- svdout$v %*% diag(svdout$d) 
#
## Run t-SNE
#source('~/Research_Local/FIt-SNE/fast_tsne.R', chdir=T)
##k <- simAlraK$k
#k <- 50
#init <- 0.0001*(obsPCA[,1:2]/sd(obsPCA[,1]))
#fitsneout <- fftRtsne( obsPCA[,1:k],initialization = init , rand_seed=3)
#df=data.frame(fitsneout, Y= as.factor(cellType))
#df <- df[!is.na(df$Y),] # Remove cells that don't have an assigned celltype
##df <- dplyr::filter(df, !df$Y %in%   c('B.Cell', 'NK.Cell', 'Neutrophil', 'Dendritic.Cell'))
#dim(df)
#(g <- ggplot(df) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+  labs(x="t-SNE 1", y="t-SNE 2"))
#ggsave(g, file="figs/test_tsne3.pdf")
#
#
#### alra_ImmGen_comparison.R
#### @Jun Zhao
#### script to generate results for method comparison on ImmGen simulation data
#
#library(ggplot2)
#library(reshape2)
#library(cowplot)
#library(rsvd)
#library(Rmagic)
#library(SAVER)
#library(scImpute)
#source("../ALRA/alra.R")
#source("alra_ImmGen_funs.R")
#
#
###=======================================================##
### Generate simulation data
#
#
### Get p for each cell type from bulk RNA-seq
#ls()
## load ImmGen bulk RNA-seq data
#load("data/ImmGen_bulk.RData")
#
## get prob
#bulkProb <- sweep(bulkCount, MARGIN = 2, STATS = colSums(bulkCount), FUN = "/")
#
#
#################################
## Simulate 
################################
### Get cell transcript count from PBMC data
#ss <- which(!bulkType %in%   c('B.Cell', 'NK.Cell', 'Neutrophil', 'Dendritic.Cell') )
##ss <-  which(!bulkType %in% c('fofo'))
#bulkType <- bulkType[ss]
#bulkCount  <- bulkCount[,ss]
#bulkNorm  <- bulkNorm[,ss]
#
## load cell size data
#load("data/PBMC_lib_size.RData")
### Simulation
#
## simulation data init
#nCell <- 30000
#nGene <- nrow(bulkCount)
#
## subsample from PBMC library size
#set.seed(3)
#simCellSize <- 0.2*sample(PBMC_lib_size, size = nCell, replace = F)
#
## simulation data matrix
#simCount <- matrix(0, nrow = nGene, ncol = nCell)
## simulate cell type
#cellType <- sample(bulkType, size = nCell, replace = T)
#table(cellType)
#
#
## get count from multinomial distribution
#for(i in 1:nCell){
#    simCount[,i] <- rmultinom(n = 1, size = simCellSize[i], prob = bulkProb[,which(bulkType == cellType[i])])
#}
#
#mean(colSums(simCount) >= 1000)
#rownames(simCount)             <- rownames(bulkCount)
#colnames(simCount)             <- paste("Cell", seq(1:nCell), sep = "_")
#
#simCount <- simCount[rowSums(simCount) > 0,]
#mean(simCount == 0) # ~0.936
#
#
## normalize data
#realCellIdx <- sample(which(colSums(simCount) >= 1), size = 5000, replace = F)
#cellType_  <- cellType[realCellIdx]
#
#length(realCellIdx)
#mean(simCount[,realCellIdx] == 0) # ~0.931
#simNorm <- log(normbycol(simCount[rowSums(simCount[,realCellIdx]) != 0,realCellIdx], 10000) + 1)
#
#dim(simCount)
#
#################################
## Look at it 
#################################
##tt <- which(!cellType_ %in%   c('B.Cell', 'NK.Cell', 'Neutrophil', 'Dendritic.Cell') )
#tt <- which(!cellType_ %in%   c('fofo') )
#length(tt)
#length(cellType_)
#dim(simNorm)
#simNorm_tt  <- simNorm[,tt]
#cellType_tt  <- cellType_[tt]
#
#row_gene_means <- rowMeans(simNorm_tt)
#simNorm_tt_c <- sweep(simNorm_tt, 1, row_gene_means, '-')
#system.time(svdout <- rsvd(simNorm_tt_c, k=100))
#obsPCA <- svdout$v %*% diag(svdout$d) 
#
## Run t-SNE
#source('~/Research_Local/FIt-SNE/fast_tsne.R', chdir=T)
##k <- simAlraK$k
#k <- 50
#init <- 0.0001*(obsPCA[,1:2]/sd(obsPCA[,1]))
#fitsneout <- fftRtsne( obsPCA[,1:k],initialization = init , rand_seed=3)
#df=data.frame(fitsneout, Y= as.factor(cellType_tt))
#df <- df[!is.na(df$Y),] # Remove cells that don't have an assigned celltype
##df <- dplyr::filter(df, !df$Y %in%   c('B.Cell', 'NK.Cell', 'Neutrophil', 'Dendritic.Cell'))
#dim(df)
#(g <- ggplot(df) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+  labs(x="t-SNE 1", y="t-SNE 2"))
#ggsave(g, file="figs/test_tsne3.pdf")
#
#
#
###=========================================================##
### Imputation
#
#
### ALRA
#set.seed(3)
##simAlraK    <- choose_k(t(as.matrix(simNorm)), q=10)
#simAlraK    <- choose_k(t(as.matrix(simNorm)))
#simAlra     <- alra(t(as.matrix(simNorm)), k = simAlraK$k)
#simNormAlra <- t(simAlra$A_norm_rank_k_cor_sc)
#
#
#
### MAGIC
#
#if(file.exists("data/simMagic.RData")){
#  load("data/simMagic.RData")
#} else {
#  reticulate::use_condaenv('magic')
#  simMagic     <- magic(data = t(as.matrix(simNorm)))
#  simNormMagic <- t(simMagic$result)
#  save(list = c("simMagic","simNormMagic"), file = "data/simMagic.RData")
#}
#
### SAVER
#
#if(file.exists("data/simSaver.RData")){
#  load("data/simSaver.RData")
#} else {
#  simSaver <- saver(x = simCount[rownames(simNorm),realCellIdx])
#  simNormSaver <- sample.saver(simSaver)
#  save(list = c("simSaver","simNormSaver"), file = "data/simSaver.RData")
#}
#simNormSaver <- log(simNormSaver + 1)
#
#
#
### scImpute
#
#if(!file.exists("data/scimpute_simulation/scimpute_count.csv")){
#    # write simulated data in file
#    print("Running scimpute")
#    write.table(
#      simCount[rowSums(simCount[,realCellIdx]) != 0,realCellIdx], "data/simulation_5000.csv", 
#      sep = ",", quote = F, row.names = T, col.names = T
#    )
#    set.seed(3)
#    scimpute_result <- scimpute(
#      count_path = "data/simulation_5000.csv",
#      out_dir = "data/scimpute_simulation/", 
#      Kcluster = 9 
#    )
#}
#scimpute_norm <- as.matrix(read.table(
#  "data/scimpute_simulation/scimpute_count.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F
#))
#scimpute_norm <- log(normbycol(scimpute_norm, 10000) + 1)
#
#
### DCA
#
#if(!file.exists("data/dca_simulation/mean.tsv")){
#    system('dca data/simulation_5000.csv data/dca_simulation/')
#}
#dca_norm <- read.table(
#    "data/dca_simulation/mean.tsv",
#    sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
#dca_norm <- as.matrix(log(normbycol(dca_norm, 10000) + 1))
#
#
#
#
#
###=======================================================##
### Evaluation
#
#
### Zero preservation
#
#dim(bulkCount)
#bulkCount[1:5,1:5]
#length(bulkType)
#bulkType
#dim(simNorm)
#simNorm[1:5,1:5]
#cellType[realCellIdx][1:5]
## original data
#zero.eval.observed <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
#                               impute.data = simNorm, impute.label = cellType[realCellIdx])
#
#
## evaluate ALRA
#zero.eval.ALRA <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
#                           impute.data = simNormAlra, impute.label = cellType[realCellIdx])
#zero.eval.ALRA[,2] <- zero.eval.ALRA[,2] - zero.eval.observed[,2]
#
#
## evaluate MAGIC
#zero.eval.MAGIC <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
#                            impute.data = simNormMagic, impute.label = cellType[realCellIdx])
#zero.eval.MAGIC[,2] <- zero.eval.MAGIC[,2] - zero.eval.observed[,2]
#
#
## evaluate SAVER
#zero.eval.SAVER <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
#                            impute.data = simNormSaver, impute.label = cellType[realCellIdx])
#zero.eval.SAVER[,2] <- zero.eval.SAVER[,2] - zero.eval.observed[,2]
#
#
## evaluate scImpute
#zero.eval.scImpute <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
#                               impute.data = scimpute_norm, impute.label = cellType[realCellIdx])
#zero.eval.scImpute[,2] <- zero.eval.scImpute[,2] - zero.eval.observed[,2]
#
#
## combine evaluation
#zero.eval.combine <- matrix(0, nrow = 4, ncol = 2)
#colnames(zero.eval.combine) <- c("y1","y2")
#rownames(zero.eval.combine) <- c("ALRA", "MAGIC", "SAVER", "scImpute")
#
#zero.eval.combine[1,] <- colSums(sweep(zero.eval.ALRA*100,1,
#                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
#zero.eval.combine[2,] <- colSums(sweep(zero.eval.MAGIC*100,1,
#                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
#zero.eval.combine[3,] <- colSums(sweep(zero.eval.SAVER*100,1,
#                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
#zero.eval.combine[4,] <- colSums(sweep(zero.eval.scImpute*100,1,
#                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
#
#zero.eval.combine.df <- data.frame(method = rownames(zero.eval.combine), zero.eval.combine, stringsAsFactors = F)
#zero.eval.combine.df[,1] <- factor(zero.eval.combine.df[,1], levels = c("MAGIC","SAVER","ALRA","scImpute"))
#
#zero.eval.combine.df$method <- factor(zero.eval.combine.df$method, levels=c("scImpute",  "SAVER", "MAGIC", "ALRA"))
#
#
## plot combined result
#gg1 <- ggplot() + theme_cowplot() + 
#  geom_bar(data = zero.eval.combine.df, stat="identity", aes(x = method, y = y1, fill = method)) + coord_flip() + 
#  geom_text(data = zero.eval.combine.df,
#            aes(x = method, y = y1, label = format(y1,digits = 2)), stat='identity', hjust = 1.3) + 
#  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
#        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,0,1,1),"cm")) + 
#  scale_y_reverse(name = "Biological zeros preserved (%)", limits = c(110,0), breaks = c(0,25,50,75,100)) + 
#  geom_hline(yintercept = 100, linetype = 2) + scale_fill_discrete(breaks=rev(levels(zero.eval.combine.df$method))) + theme(legend.position = "none")
#gg2 <- ggplot() + theme_cowplot() + 
#  geom_bar(data = zero.eval.combine.df, stat="identity", aes(x = method, y = y2, fill = method)) + coord_flip() + 
#  geom_text(data = zero.eval.combine.df,
#            aes(x = method, y = y2, label = format(y2,digits = 2)), stat='identity', hjust = -0.3) + 
#  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
#        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,1,1,0),"cm")) + 
#  scale_y_continuous(name = "Technical zeros completed (%)", limits = c(0,110), breaks = c(0,25,50,75,100)) + 
#  geom_hline(yintercept = 100, linetype = 2) + theme(legend.position = "none")
#(g <- gridExtra::grid.arrange(gg1,gg2,nrow=1))
#ggsave(g, filename = "figs/ImmGen_zero_preserve_4method_10000.pdf", width = 10, height = 2.5)
#
## legend
#gg1 <- gg1 + theme(legend.position = "bottom", legend.title = element_blank()) 
#ggsave(grid::grid.draw(g_legend(gg1)), filename = "figs/ImmGen_zero_preserve_4method_legend.pdf", width = 4, height = 1)
#
#zero.eval.combinne
#
#
### Correlation with bulk RNA-seq
#
## observed
#cor.eval.observed <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = as.matrix(simNorm),
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.observed <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = as.matrix(simNorm),
#    impute.label = cellType[realCellIdx]
#)
#
#
## ALRA
#cor.eval.ALRA <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormAlra,
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.ALRA <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormAlra,
#    impute.label = cellType[realCellIdx]
#)
#
#
## MAGIC
#cor.eval.MAGIC <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormMagic,
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.MAGIC <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormMagic,
#    impute.label = cellType[realCellIdx]
#)
#
#
## SAVER
#cor.eval.SAVER <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormSaver,
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.SAVER <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = simNormSaver,
#    impute.label = cellType[realCellIdx]
#)
#
#
## scImpute
#cor.eval.scImpute <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = scimpute_norm,
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.scImpute <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = scimpute_norm,
#    impute.label = cellType[realCellIdx]
#)
#
#
## DCA
#cor.eval.DCA <- corEval(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = dca_norm,
#    impute.label = cellType[realCellIdx]
#)
#cor.cell.DCA <- corPerCell(
#    true.data = bulkNorm[rownames(simNorm),],
#    true.label = bulkType,
#    impute.data = dca_norm,
#    impute.label = cellType[realCellIdx]
#)
#
## plot 4 correlation diff together
#gg1 <- plotCellCorDiff(
#    cor.cell.observed,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + annotate(
#    "text",
#    x = 0.75,
#    y = -0.08,
#    label = paste("Error: ", format(
#        (1-sum(diag(cor.eval.observed))/sum(cor.eval.observed))*100,
#        nsmall = 2
#    ), "%", sep = "")
#)
#gg2 <- plotCellCorDiff(
#    cor.cell.ALRA,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + 
#  annotate(
#      "text",
#      x = 0.75,
#      y = -0.08,
#      label = paste(
#          "Error: ",
#          format((1-sum(diag(cor.eval.ALRA))/sum(cor.eval.ALRA))*100,nsmall = 2),
#          "%",
#          sep = ""
#      )
#  )
#gg3 <- plotCellCorDiff(
#    cor.cell.MAGIC,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + annotate(
#    "text",
#    x = 0.72,
#    y = -0.08,
#    label = paste(
#        "Error: ",
#        format((1-sum(diag(cor.eval.MAGIC))/sum(cor.eval.MAGIC))*100,nsmall = 2),
#        "%",
#        sep = ""
#    )
#)
#gg4 <- plotCellCorDiff(
#    cor.cell.SAVER,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + annotate(
#    "text",
#    x = 0.72,
#    y = -0.08,
#    label = paste(
#        "Error: ",
#        format((1-sum(diag(cor.eval.SAVER))/sum(cor.eval.SAVER))*100,nsmall = 2),
#        "%",
#        sep = ""
#    )
#)
#gg5 <- plotCellCorDiff(
#    cor.cell.scImpute,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + annotate(
#    "text",
#    x = 0.72,
#    y = -0.08,
#    label = paste(
#        "Error: ",
#        format((1-sum(diag(cor.eval.scImpute))/sum(cor.eval.scImpute))*100,nsmall = 2),
#        "%",
#        sep = ""
#    )
#)
#gg6 <- plotCellCorDiff(
#    cor.cell.DCA,
#    title = "",
#    input.cell.label = cellType[realCellIdx]
#) + 
#  annotate(
#      "text",
#      x = 0.72,
#      y = -0.08,
#      label = paste("Error: ", format((1-sum(diag(cor.eval.DCA))/sum(cor.eval.DCA))*100,nsmall = 2), "%", sep = "")
#  )
#
#
#g <- gridExtra::grid.arrange(gg1,gg2,gg3,gg4,gg5,gg6, nrow = 1)
#ggsave(g, file='figs/ImmGen_cell_correlationDiff_5000_5method_1row.pdf', width=15, height=2.5)
#
#
## plot legend
#ggsave(grid::grid.draw(g_legend(plotCellCor(cor.cell.observed, title = "Observed", 
#                                            input.cell.label = gsub(".", " ", cellType[realCellIdx], fixed = T), 
#                                            legend = T))), 
#       filename = "figs/ImmGen_correlation_legend.pdf", width = 8, height = 2)
#
#
## plot mean in the same canvas
#cor.diff.combine <- cbind(by(cor.cell.observed[,1]-cor.cell.observed[,2], INDICES = cellType[realCellIdx], mean),
#                          by(cor.cell.ALRA[,1]-cor.cell.ALRA[,2], INDICES = cellType[realCellIdx], mean),
#                          by(cor.cell.MAGIC[,1]-cor.cell.MAGIC[,2], INDICES = cellType[realCellIdx], mean),
#                          by(cor.cell.SAVER[,1]-cor.cell.SAVER[,2], INDICES = cellType[realCellIdx], mean),
#                          by(cor.cell.scImpute[,1]-cor.cell.scImpute[,2], INDICES = cellType[realCellIdx], mean),
#                          by(cor.cell.DCA[,1]-cor.cell.DCA[,2], INDICES = cellType[realCellIdx], mean))
#colnames(cor.diff.combine) <- c("Observed","ALRA","MAGIC","SAVER","scImpute","DCA")
#
#gg1 <- ggplot(data = melt(cor.diff.combine), aes(x = Var2, y = value, group = Var1, col = Var1)) + theme_cowplot() + 
#  geom_line() + 
#  geom_point() + 
#  xlab("") + ylab("") + ggtitle("") + 
#  scale_x_discrete(labels = c("","","","","","")) + 
#  theme(legend.position = "none")
#
#ggsave(gg1, filename = "figs/ImmGen_cell_correlationDiffLine_5000_5method.pdf", height = 2.5, width = 3)
#
#
#    # original data
#    #zero.eval.observed <- zeroEval(true.data = bulkCount[keep.genes,] , true.label = colnames(bulkCount), 
#    #                               impute.data = simNorm, impute.label = cellType.sub)
#    ## evaluate ALRA
#    #zero.eval.ALRA <- zeroEval(true.data = bulkCount[keep.genes,] , true.label = colnames(bulkCount), 
#    #                            impute.data = simNormAlra, impute.label = cellType.sub)
#    ## gvaluate MAGIC
#    #zero.eval.MAGIC <- zeroEval(true.data = bulkCount[keep.genes,] , true.label = colnames(bulkCount), 
#    #                            impute.data = simNormMagic, impute.label = cellType.sub)
#    #zero.eval.MAGIC[,2] <- zero.eval.MAGIC[,2] - zero.eval.observed[,2]
#    ## combine evaluation
#    #zero.eval.combine <- matrix(0, nrow = 4, ncol = 2)
#    #colnames(zero.eval.combine) <- c("y1","y2")
#    #rownames(zero.eval.combine) <- c("ALRA", "MAGIC", "SAVER", "scImpute")
#    #alra.result <- colSums(sweep(zero.eval.ALRA*100,1, table(cellType)[rownames(zero.eval.observed)]/length(cellType),"*"))
#    #magic.result <- colSums(sweep(zero.eval.MAGIC*100,1, table(cellType)[rownames(zero.eval.observed)]/length(cellType),"*"))
#    #true.zeros[fraci,1]  <- alra.result[1]
#    #true.nonzeros[fraci,1]  <- alra.result[2]
#    #true.zeros[fraci,2]  <- magic.result[1]
#    #true.nonzeros[fraci,2]  <- magic.result[2]
