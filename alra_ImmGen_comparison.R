### alra_ImmGen_comparison.R
### @Jun Zhao
### script to generate results for method comparison on ImmGen simulation data

library(ggplot2)
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
load("data/ImmGen_bulk.RData")

# get prob
bulkProb <- sweep(bulkCount, MARGIN = 2, STATS = colSums(bulkCount), FUN = "/")



if(file.exists("data/simData.RData")){
  print("Loading...")
  load("data/simData.RData")
} else {
  ## Get cell transcript count from PBMC data
  
  # load cell size data
  load("data/PBMC_lib_size.RData")
  
  
  ## Simulation
  
  # simulation data init
  nCell <- 50000
  nGene <- nrow(bulkCount)
  
  # subsample from PBMC library size
  set.seed(3)
  simCellSize <- 0.2*sample(PBMC_lib_size, size = nCell, replace = F)
  
  # simulation data matrix
  simCount <- matrix(0, nrow = nGene, ncol = nCell)
  # simulate cell type
  cellType <- sample(bulkType, size = nCell, replace = T)
  table(cellType)
  
  
  # get count from multinomial distribution
  for(i in 1:nCell){
    simCount[,i] <- rmultinom(n = 1, size = simCellSize[i], prob = bulkProb[,which(bulkType == cellType[i])])
  }
  
  mean(colSums(simCount) >= 1000)
  rownames(simCount)             <- rownames(bulkCount)
  colnames(simCount)             <- paste("Cell", seq(1:nCell), sep = "_")
  
  simCount <- simCount[rowSums(simCount) > 0,]
  mean(simCount == 0) # ~0.936
  
  
  # normalize data
  realCellIdx <- sample(which(colSums(simCount) >= 1), size = 10000, replace = F)
  mean(simCount[,realCellIdx] == 0) # ~0.931
  simNorm <- log(normbycol(simCount[rowSums(simCount[,realCellIdx]) != 0,realCellIdx], 10000) + 1)
  save(simCount, cellType, simNorm, realCellIdx, file = "data/simData.RData")

}





##=========================================================##
## Imputation


## ALRA
set.seed(3)
#simAlraK    <- choose_k(t(as.matrix(simNorm)), q=10)
simAlraK    <- choose_k(t(as.matrix(simNorm)))
simAlra     <- alra(t(as.matrix(simNorm)), k = simAlraK$k)
simNormAlra <- t(simAlra$A_norm_rank_k_cor_sc)



## MAGIC

if(file.exists("data/simMagic.RData")){
  load("data/simMagic.RData")
} else {
  reticulate::use_condaenv('magic')
  simMagic     <- magic(data = t(as.matrix(simNorm)))
  simNormMagic <- t(simMagic$result)
  save(list = c("simMagic","simNormMagic"), file = "data/simMagic.RData")
}

## SAVER

if(file.exists("data/simSaver.RData")){
  load("data/simSaver.RData")
} else {
  simSaver <- saver(x = simCount[rownames(simNorm),realCellIdx])
  simNormSaver <- sample.saver(simSaver)
  save(list = c("simSaver","simNormSaver"), file = "data/simSaver.RData")
}
simNormSaver <- log(simNormSaver + 1)



## scImpute

if(!file.exists("data/scimpute_simulation/scimpute_count.csv")){
    # write simulated data in file
    print("Running scimpute")
    write.table(
      simCount[rowSums(simCount[,realCellIdx]) != 0,realCellIdx], "data/simulation_5000.csv", 
      sep = ",", quote = F, row.names = T, col.names = T
    )
    set.seed(3)
    scimpute_result <- scimpute(
      count_path = "data/simulation_5000.csv",
      out_dir = "data/scimpute_simulation/", 
      Kcluster = 9 
    )
}
scimpute_norm <- as.matrix(read.table(
  "data/scimpute_simulation/scimpute_count.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F
))
scimpute_norm <- log(normbycol(scimpute_norm, 10000) + 1)


## DCA

if(!file.exists("data/dca_simulation/mean.tsv")){
    system('dca data/simulation_5000.csv data/dca_simulation/')
}
dca_norm <- read.table(
    "data/dca_simulation/mean.tsv",
    sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
dca_norm <- as.matrix(log(normbycol(dca_norm, 10000) + 1))





##=======================================================##
## Evaluation


## Zero preservation

# original data
zero.eval.observed <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
                               impute.data = simNorm, impute.label = cellType[realCellIdx])


# evaluate ALRA
zero.eval.ALRA <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
                           impute.data = simNormAlra, impute.label = cellType[realCellIdx])
zero.eval.ALRA[,2] <- zero.eval.ALRA[,2] - zero.eval.observed[,2]


# evaluate MAGIC
zero.eval.MAGIC <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
                            impute.data = simNormMagic, impute.label = cellType[realCellIdx])
zero.eval.MAGIC[,2] <- zero.eval.MAGIC[,2] - zero.eval.observed[,2]


# evaluate SAVER
zero.eval.SAVER <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
                            impute.data = simNormSaver, impute.label = cellType[realCellIdx])
zero.eval.SAVER[,2] <- zero.eval.SAVER[,2] - zero.eval.observed[,2]


# evaluate scImpute
zero.eval.scImpute <- zeroEval(true.data = bulkCount[rownames(simNorm),], true.label = bulkType, 
                               impute.data = scimpute_norm, impute.label = cellType[realCellIdx])
zero.eval.scImpute[,2] <- zero.eval.scImpute[,2] - zero.eval.observed[,2]


# combine evaluation
zero.eval.combine <- matrix(0, nrow = 4, ncol = 2)
colnames(zero.eval.combine) <- c("y1","y2")
rownames(zero.eval.combine) <- c("ALRA", "MAGIC", "SAVER", "scImpute")

zero.eval.combine[1,] <- colSums(sweep(zero.eval.ALRA*100,1,
                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
zero.eval.combine[2,] <- colSums(sweep(zero.eval.MAGIC*100,1,
                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
zero.eval.combine[3,] <- colSums(sweep(zero.eval.SAVER*100,1,
                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))
zero.eval.combine[4,] <- colSums(sweep(zero.eval.scImpute*100,1,
                                       table(cellType[realCellIdx])[rownames(zero.eval.observed)]/length(realCellIdx),"*"))

zero.eval.combine.df <- data.frame(method = rownames(zero.eval.combine), zero.eval.combine, stringsAsFactors = F)
zero.eval.combine.df[,1] <- factor(zero.eval.combine.df[,1], levels = c("MAGIC","SAVER","ALRA","scImpute"))

zero.eval.combine.df$method <- factor(zero.eval.combine.df$method, levels=c("scImpute",  "SAVER", "MAGIC", "ALRA"))


# plot combined result
gg1 <- ggplot() + theme_cowplot() + 
  geom_bar(data = zero.eval.combine.df, stat="identity", aes(x = method, y = y1, fill = method)) + coord_flip() + 
  geom_text(data = zero.eval.combine.df,
            aes(x = method, y = y1, label = format(y1,digits = 2)), stat='identity', hjust = 1.3) + 
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,0,1,1),"cm")) + 
  scale_y_reverse(name = "Biological zeros preserved (%)", limits = c(110,0), breaks = c(0,25,50,75,100)) + 
  geom_hline(yintercept = 100, linetype = 2) + scale_fill_discrete(breaks=rev(levels(zero.eval.combine.df$method))) + theme(legend.position = "none")
gg2 <- ggplot() + theme_cowplot() + 
  geom_bar(data = zero.eval.combine.df, stat="identity", aes(x = method, y = y2, fill = method)) + coord_flip() + 
  geom_text(data = zero.eval.combine.df,
            aes(x = method, y = y2, label = format(y2,digits = 2)), stat='identity', hjust = -0.3) + 
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,1,1,0),"cm")) + 
  scale_y_continuous(name = "Technical zeros completed (%)", limits = c(0,110), breaks = c(0,25,50,75,100)) + 
  geom_hline(yintercept = 100, linetype = 2) + theme(legend.position = "none")
(g <- gridExtra::grid.arrange(gg1,gg2,nrow=1))
ggsave(g, filename = "figs/ImmGen_zero_preserve_4method_10000.pdf", width = 10, height = 2.5)

# legend
gg1 <- gg1 + theme(legend.position = "bottom", legend.title = element_blank()) 
ggsave(grid::grid.draw(g_legend(gg1)), filename = "figs/ImmGen_zero_preserve_4method_legend.pdf", width = 4, height = 1)

zero.eval.combinne


## Correlation with bulk RNA-seq

# observed
cor.eval.observed <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = as.matrix(simNorm),
    impute.label = cellType[realCellIdx]
)
cor.cell.observed <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = as.matrix(simNorm),
    impute.label = cellType[realCellIdx]
)


# ALRA
cor.eval.ALRA <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormAlra,
    impute.label = cellType[realCellIdx]
)
cor.cell.ALRA <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormAlra,
    impute.label = cellType[realCellIdx]
)


# MAGIC
cor.eval.MAGIC <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormMagic,
    impute.label = cellType[realCellIdx]
)
cor.cell.MAGIC <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormMagic,
    impute.label = cellType[realCellIdx]
)


# SAVER
cor.eval.SAVER <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormSaver,
    impute.label = cellType[realCellIdx]
)
cor.cell.SAVER <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = simNormSaver,
    impute.label = cellType[realCellIdx]
)


# scImpute
cor.eval.scImpute <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = scimpute_norm,
    impute.label = cellType[realCellIdx]
)
cor.cell.scImpute <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = scimpute_norm,
    impute.label = cellType[realCellIdx]
)


# DCA
cor.eval.DCA <- corEval(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = dca_norm,
    impute.label = cellType[realCellIdx]
)
cor.cell.DCA <- corPerCell(
    true.data = bulkNorm[rownames(simNorm),],
    true.label = bulkType,
    impute.data = dca_norm,
    impute.label = cellType[realCellIdx]
)

# plot 4 correlation diff together
gg1 <- plotCellCorDiff(
    cor.cell.observed,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + annotate(
    "text",
    x = 0.75,
    y = -0.08,
    label = paste("Error: ", format(
        (1-sum(diag(cor.eval.observed))/sum(cor.eval.observed))*100,
        nsmall = 2
    ), "%", sep = "")
)
gg2 <- plotCellCorDiff(
    cor.cell.ALRA,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + 
  annotate(
      "text",
      x = 0.75,
      y = -0.08,
      label = paste(
          "Error: ",
          format((1-sum(diag(cor.eval.ALRA))/sum(cor.eval.ALRA))*100,nsmall = 2),
          "%",
          sep = ""
      )
  )
gg3 <- plotCellCorDiff(
    cor.cell.MAGIC,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + annotate(
    "text",
    x = 0.72,
    y = -0.08,
    label = paste(
        "Error: ",
        format((1-sum(diag(cor.eval.MAGIC))/sum(cor.eval.MAGIC))*100,nsmall = 2),
        "%",
        sep = ""
    )
)
gg4 <- plotCellCorDiff(
    cor.cell.SAVER,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + annotate(
    "text",
    x = 0.72,
    y = -0.08,
    label = paste(
        "Error: ",
        format((1-sum(diag(cor.eval.SAVER))/sum(cor.eval.SAVER))*100,nsmall = 2),
        "%",
        sep = ""
    )
)
gg5 <- plotCellCorDiff(
    cor.cell.scImpute,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + annotate(
    "text",
    x = 0.72,
    y = -0.08,
    label = paste(
        "Error: ",
        format((1-sum(diag(cor.eval.scImpute))/sum(cor.eval.scImpute))*100,nsmall = 2),
        "%",
        sep = ""
    )
)
gg6 <- plotCellCorDiff(
    cor.cell.DCA,
    title = "",
    input.cell.label = cellType[realCellIdx]
) + 
  annotate(
      "text",
      x = 0.72,
      y = -0.08,
      label = paste("Error: ", format((1-sum(diag(cor.eval.DCA))/sum(cor.eval.DCA))*100,nsmall = 2), "%", sep = "")
  )


g <- gridExtra::grid.arrange(gg1,gg2,gg3,gg4,gg5,gg6, nrow = 1)
ggsave(g, file='figs/ImmGen_cell_correlationDiff_5000_5method_1row.pdf', width=15, height=2.5)


# plot legend
ggsave(grid::grid.draw(g_legend(plotCellCor(cor.cell.observed, title = "Observed", 
                                            input.cell.label = gsub(".", " ", cellType[realCellIdx], fixed = T), 
                                            legend = T))), 
       filename = "figs/ImmGen_correlation_legend.pdf", width = 8, height = 2)


# plot mean in the same canvas
cor.diff.combine <- cbind(by(cor.cell.observed[,1]-cor.cell.observed[,2], INDICES = cellType[realCellIdx], mean),
                          by(cor.cell.ALRA[,1]-cor.cell.ALRA[,2], INDICES = cellType[realCellIdx], mean),
                          by(cor.cell.MAGIC[,1]-cor.cell.MAGIC[,2], INDICES = cellType[realCellIdx], mean),
                          by(cor.cell.SAVER[,1]-cor.cell.SAVER[,2], INDICES = cellType[realCellIdx], mean),
                          by(cor.cell.scImpute[,1]-cor.cell.scImpute[,2], INDICES = cellType[realCellIdx], mean),
                          by(cor.cell.DCA[,1]-cor.cell.DCA[,2], INDICES = cellType[realCellIdx], mean))
colnames(cor.diff.combine) <- c("Observed","ALRA","MAGIC","SAVER","scImpute","DCA")

gg1 <- ggplot(data = melt(cor.diff.combine), aes(x = Var2, y = value, group = Var1, col = Var1)) + theme_cowplot() + 
  geom_line() + 
  geom_point() + 
  xlab("") + ylab("") + ggtitle("") + 
  scale_x_discrete(labels = c("","","","","","")) + 
  theme(legend.position = "none")

ggsave(gg1, filename = "figs/ImmGen_cell_correlationDiffLine_5000_5method.pdf", height = 2.5, width = 3)


