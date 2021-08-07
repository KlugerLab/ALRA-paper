# The down-sample analysis was inspired by Huang et al.
# (2018). In order to compare with their results, we closely followed their
# analysis, the code for which is available at
# https://github.com/mohuangx/SAVER-paper, and data available at
# https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0


#setwd("/data/jz437/Drop/ALRA/downsample/")
#save(list = ls(), file = "./alra_downsample.RData")
library(Seurat)
library(cowplot)
library(ranger)
library(ggplot2)

source('convenience.R')
set.seed(3)
source("../ALRA/alra.R")


## Prepare data

# please download data with "ref" "samp.rds" "samp_saver.rds" "samp_magic.rds" "samp_scimpute.rds" patterns

# load data
ref.files <- list.files("data", pattern = "ref", full.names = TRUE)
samp.files <- list.files("data", pattern = "samp.rds", full.names = TRUE)
saver.files <- list.files("data", pattern = "samp_saver.rds", full.names = TRUE)
magic.files <- list.files("data", pattern = "samp_magic.rds", full.names = TRUE)
scimpute.files <- list.files("data", pattern = "samp_scimpute.rds", full.names = TRUE)

dat <- vector("list", 4)
names(dat) <- c("baron", "chen", "manno", "zeisel")

normalizeData <- function(x, y = x) {
  sf <- colSums(y)/mean(colSums(y))
  return(sweep(x, 2, sf, "/"))
}

for (i in 1:4) {
  dat[[i]] <- vector("list", 9)
  names(dat[[i]]) <- c("X", "Y", "saver", "alra", "dca", "X.norm", "Y.norm", "magic", "scimpute")
  dat[[i]][[1]] <- readRDS(ref.files[i])
  dat[[i]][[2]] <- readRDS(samp.files[i])
  dat[[i]][[3]] <- readRDS(saver.files[i])
  dat[[i]][[8]] <- readRDS(magic.files[i])
  dat[[i]][[9]] <- readRDS(scimpute.files[i])
  dat[[i]][[9]][is.na(dat[[i]][[9]])] <- 0
}

n.cells <- sapply(dat, function(x) ncol(x[[1]]))
n.genes <- sapply(dat, function(x) nrow(x[[1]]))
cell.names <- sapply(dat, function(x) colnames(x[[1]]))
gene.names <- sapply(dat, function(x) rownames(x[[1]]))

colnames(dat[[1]][[9]]) <- gsub(".","-", colnames(dat[[1]][[9]]),fixed=T)
colnames(dat[[3]][[9]]) <- substring(colnames(dat[[3]][[9]]),2)
colnames(dat[[4]][[9]]) <- substring(colnames(dat[[4]][[9]]),2)
for (i in 1:4) {
  dat[[i]][[6]] <- normalizeData(dat[[i]][[1]])
  dat[[i]][[7]] <- normalizeData(dat[[i]][[2]])
  dat[[i]][[9]] <- normalizeData(dat[[i]][[9]])
}
str(dat)


# alra
alraResult <- function(x){
  x.norm <- normalize_data(t(x))
  x.alra <- t(alra(x.norm)[[3]])
  rownames(x.alra) <- rownames(x)
  colnames(x.alra) <- colnames(x)
  return(x.alra)
}

for(i in 1:4){
  dat[[i]][[4]] <- alraResult(dat[[i]][[2]])
}


## choose k

 #### calculate PCs using jackstraw
 #calc.jackstraw <- function(x, jack = TRUE) {
 #  x.seurat <- CreateSeuratObject(raw.data = x)
 #  x.seurat <- NormalizeData(x.seurat)
 #  x.seurat <- ScaleData(x.seurat)
 #  x.seurat <- FindVariableGenes(x.seurat, do.plot = FALSE)
 #  x.seurat <- RunPCA(x.seurat, pc.genes = x.seurat@var.genes,
 #                     pcs.compute = 25, do.print = FALSE)
 #  if (jack) {
 #    x.seurat <- JackStraw(x.seurat, num.pc = 25, do.print = TRUE)
 #  }
 #  return(x.seurat)
 #}
 #
 #dat.seurat <- vector("list", 4)
 #names(dat.seurat) <- names(dat)
 #
 #for (i in 1:4) {
 #  dat.seurat[[i]] <- vector("list", 9)
 #  names(dat.seurat[[i]]) <- c("X", "Y", "saver", "alra", "dca", "X.norm", "Y.norm", "magic", "scimpute")
 #  for (j in c("X", "Y", "saver", "alra", "dca", "magic", "scimpute")) {
 #    dat.seurat[[i]][[j]] <- calc.jackstraw(dat[[i]][[j]], jack = FALSE)
 #    print(j)
 #  }
 #}
 #
 #
 #for (i in 1:4) {
 #  for (j in 1:5) {
 #    dat.seurat[[i]][[j]]@raw.data <- NULL
 #  }
 #}
 #
 #jackstraw.results <- vector("list", 4)
 #
 #for (i in 1:4) {
 #  jackstraw.results[[i]] <- sapply(dat.seurat[[i]], function(x) 
 #    levels(JackStrawPlot(x, PCs = 1:25)$data$PC.Score))
 #}
 




#alraChooseK <- function(x){
#  x.norm <- normalize_data(t(x))
#  return(choose_k(x.norm))
#}
#dat.k <- vector("list",4)
#names(dat.k) <- names(dat)
#for(i in 1:4){
#  dat.k[[i]] <- alraChooseK(dat[[i]][[2]])
#}
#
#dat.k[[2]]$k
#plot(dat.k[[2]]$d[4:100])
#plot(diff(dat.k[[2]]$d)[4:100])
#plot(dat.k[[2]]$pvals)



# DCA
dcaResult <- function(x, name, do.read = F){
  if(do.read){
    x.dca <- read.table(
      paste("data/DCA_", name, "/dca_default/mean_norm.tsv", sep = ""), sep="\t", header = TRUE, row.names = 1
    )
    return(as.matrix(x.dca))
  }
  dir.create(paste("data/DCA_", name, sep = ""))
  write.csv(x, file = paste("data/DCA_", name, "/count.csv", sep = ""), row.names = T, col.names = T, quote = F)
  command <- paste("dca data/DCA_", name, "/count.csv data/DCA_", name, "/dca_default", sep = "")
  print(command)
  system(command)
  return(NULL)
}

for(i in 1:4){
  dcaResult(dat[[i]][[2]], name = names(dat)[i])
}
for(i in 1:4){
  dat[[i]][[5]] <- dcaResult(dat[[i]][[2]], name = names(dat)[i], do.read = T)
}



fofo


###############################################################################
## Cell clustering


## Get t-sne embedding

# function with normalized data
get.tsne <- function(X, n.pc, do.norm = F){
  X.S <- CreateSeuratObject(raw.data = X)
  X.S@data <- X
  if(do.norm){
    X.S <- NormalizeData(X.S)
  }
  X.S <- ScaleData(X.S)
  X.S <- FindVariableGenes(X.S, do.plot = F)
  X.S <- RunPCA(X.S, pcs.compute = n.pc, do.print = F)
  X.S <- RunTSNE(X.S, dims.use = 1:n.pc, check_duplicates = F, seed.use = 3)
  return(X.S@dr$tsne@cell.embeddings)
}


# function to plot tsne embedding
plot.tsne <- function(x, label, size = 0.5, cols.use = NULL, title = NULL, titlesize=10){
  myggplot <- ggplot() + theme_cowplot() + 
    geom_point(data = data.frame(x = x[,1], y = x[,2], class = as.factor(label)), aes(x, y, col = class), size = size) + theme_cowplot() +
    theme(axis.title = element_blank(), axis.text = element_blank(), 
                  axis.ticks = element_blank(), 
                          plot.title = element_text(size = titlesize, face = "plain"), 
                          legend.position="none", 
                          plot.margin=unit(c(0.05,0.05,0.05,0.05 ), "cm"),
                        axis.line.y.left=element_line(colour="black",size=0.2),
                    axis.line.x.bottom=element_line(colour="black",size=0.2))
    #theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank()) + theme_cowplot()
  if(!is.null(cols.use)){
    myggplot <- myggplot + 
      scale_color_manual(values = cols.use)
  }
  if(!is.null(title)){
    myggplot <- myggplot + ggtitle(title)
  }
  return(myggplot)
}


# load info
info <- readRDS("data/fig2d_tsne.rds")


# tsne info for each method
tsne <- vector("list", length = 4)
names(tsne) <- names(dat)
for(i in 1:4){
  tsne[[i]] <- vector("list", length = 6)
  names(tsne[[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
}

for(i in 1:4){
  for(j in 1:6){
    ind <- c(7,3,4,5,8,9)
    if(j == 3 | j == 5){
      tsne[[i]][[j]] <- get.tsne(X = dat[[i]][[ind[j]]], n.pc = dat.k[[i]]$k)
    } else if (j == 2) {
      tsne[[i]][[j]] <- get.tsne(X = log(dat[[i]][[ind[j]]]$estimate + 1), n.pc = dat.k[[i]]$k)
    } else {
      tsne[[i]][[j]] <- get.tsne(X = log(dat[[i]][[ind[j]]] + 1), n.pc = dat.k[[i]]$k)
    }
  }
}
for(i in 1:4){
  rownames(tsne[[i]][[4]]) <- colnames(dat[[i]][[1]])
}



## Get OOB error

# function to get var genes in Seurat
get.vargenes <- function(X, do.norm = F){
  X.S <- CreateSeuratObject(raw.data = X)
  X.S@data <- X
  if(do.norm){
    X.S <- NormalizeData(X.S)
  }
  X.S <- ScaleData(X.S)
  X.S <- FindVariableGenes(X.S, do.plot = F)
  return(X.S@var.genes)
}


# var genes for each method
vargenes <- vector("list", length = 4)
names(vargenes) <- names(dat)
for(i in 1:4){
  vargenes[[i]] <- vector("list", length = 7)
  names(vargenes[[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute", "Ref")
}

for(i in 1:4){
  for(j in 1:7){
    ind <- c(7,3,4,5,8,9, 1)
    if(j == 3 | j == 5){
      vargenes[[i]][[j]] <- get.vargenes(X = dat[[i]][[ind[j]]])
    } else if (j == 2) {
      vargenes[[i]][[j]] <- get.vargenes(X = log(dat[[i]][[ind[j]]]$estimate + 1))
    } else {
      vargenes[[i]][[j]] <- get.vargenes(X = log(dat[[i]][[ind[j]]] + 1))
    }
  }
}


# function to get OOB error
get.rfr <- function(exprs.mat, Y){
  df <- data.frame(t(exprs.mat), Y = as.factor(Y))
  rfr <- ranger(Y ~ ., df)
  return(rfr)
}
# get OOB error for each dataset and each method
oob.error <- vector("list", length = 4)
names(oob.error) <- names(dat)

for(i in 1:4){
  oob.error[[i]] <- vector("list", length = 7)
  names(oob.error[[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute", "Ref")
}
for(i in 1:4){
  for(j in 1:7){
    ind <- c(7,3,4,5,8,9,1)
    if(j == 3 | j == 5){
      dat.j <- dat[[i]][[ind[j]]]
    } else if (j == 2) {
      dat.j <- log(dat[[i]][[ind[j]]]$estimate + 1)
    } else {
      dat.j <- log(dat[[i]][[ind[j]]] + 1)
    }
    
    if(j == 4 | j == 6){
      colnames(dat.j) <- colnames(dat[[i]][[1]])
    }
    oob.error[[i]][[j]] <- get.rfr(
      exprs.mat = dat.j[vargenes[[i]][[j]],names(info[[2]][[i]]$Ref)], Y = info[[2]][[i]]$Ref
    )
  }
}


## Plot

# t-SNE plot
tsne.cols <- list(
  c("#8dd3c7", "#ffd92f", "#e78ac3", "#bebada", "#80b1d3", "#fc8d62", "#b3de69"), 
  NULL, 
  c("#8dd3c7", "#ffd92f", "#e78ac3", "#bebada", "#80b1d3", "#fc8d62", "#b3de69"), 
  c("#8dd3c7", "#ffd92f", "#e78ac3", "#bebada", "#80b1d3", "#fc8d62", "#b3de69")
)

for(i in 1:4){
  ggs <- list()
  ggs[[1]] <- plot.tsne(
    x = info[[1]][[i]]$Ref, label = info[[2]][[i]]$Ref, 
    cols.use = tsne.cols[[i]], title = "Ref", size=0.1
  ) + annotation_compass(sprintf("Error: %.1f%%", oob.error[[i]][[7]]$prediction.error*100), "SE",fontsize=9 )
  jj <- 2
  for(j in c(1,3,5,2,6,4)){
    ggs[[jj]] <- plot.tsne(
      x = tsne[[i]][[j]][names(info[[2]][[i]]$Ref),], label = info[[2]][[i]]$Ref, 
      size=0.1,
      cols.use = tsne.cols[[i]], title = names(tsne[[1]])[j]
    ) +
  annotation_compass(sprintf("Error: %.1f%%", oob.error[[i]][[j]]$prediction.error*100), "SE",fontsize=9 )
    jj <- jj + 1
  }
  g <- gridExtra::arrangeGrob(grobs = ggs[-1], nrow = 2)
  g_<- gridExtra::grid.arrange(ggs[[1]], g,
               widths=c(1,3),
                layout_matrix=  rbind(c(NA, 2),
                                      c(1, 2),
                                      c(1,2),
                                      c(NA,2)))
    ggsave(plot=g_,filename=sprintf('figs/%s_tsne_5method.pdf', names(dat)[[i]]),width=5.8,height=2.9)
}

#gridExtra::grid.arrange(grobs = ggs, ncol = 2)




## In the following sections, correlation analysis is done, 
## results not shown in the paper, maybe we should delete these?

###############################################################################
## Correlation with reference plots

# calculate correlation
get.cor.gene <- function(X, Y, method = "spearman") {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ], method = method))
}

get.cor.cell <- function(X, Y, method = "spearman") {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i], method = method))
}

cor.dat <- vector("list", 2)
names(cor.dat) <- c("gene", "cell")

for (i in 1:2) {
  cor.dat[[i]] <- vector("list", 4)
  names(cor.dat[[i]]) <- names(dat)
  for (j in 1:4) {
    cor.dat[[i]][[j]] <- vector("list", 4)
    names(cor.dat[[i]][[j]]) <- c("Obs","SAVER","ALRA","DCA")
  }
}

for (i in 1:4) {
  for (j in 1:4) {
    ind <- c(7,3,4,5)
    if(j != 2){
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]])
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]])
    } else {
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
    }
  }
}

cor.pearson.dat <- cor.dat
for (i in 1:4) {
  for (j in 1:4) {
    ind <- c(7,3,4,5)
    if(j != 2){
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]], method = "pearson")
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]], method = "pearson")
    } else {
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate, method = "pearson")
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate, method = "pearson")
    }
  }
}


# scatter plot with cor
plot.cor <- function(x, y, xlab = "", ylab = "", method = "pearson"){
  mycor <- cor(x,y, method = method)
  myggplot <- ggplot() + theme_cowplot() + 
    geom_point(
      data = data.frame(x = x, y = y),
      aes(x,y)
    ) + xlab(xlab) + ylab(ylab) + ggtitle(format(mycor, digits = 5))
  
  return(myggplot)
}


# cell-cell correlation scatter plot

datoi <- 2
coi <- 4038

# pearson
g1 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[7]][,coi],
  xlab = "Ref", ylab = "Obs"
)

g2 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[7]][,coi] + 1),
  xlab = "Ref (log)", ylab = "Obs (log)"
)

g3 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[3]]$estimate[,coi],
  xlab = "Ref", ylab = "SAVER"
)

g4 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[3]]$estimate[,coi] + 1),
  xlab = "Ref (log)", ylab = "SAVER (log)"
)

g5 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[4]][,coi],
  xlab = "Ref", ylab = "ALRA"
)

g6 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = dat[[datoi]][[4]][,coi],
  xlab = "Ref (log)", ylab = "ALRA"
)

g7 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[5]][,coi],
  xlab = "Ref", ylab = "DCA"
)

g8 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[5]][,coi]),
  xlab = "Ref (log)", ylab = "DCA (log)"
)

gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, ncol = 4)

# spearman
g1 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[7]][,coi],
  xlab = "Ref", ylab = "Obs", method = "spearman"
)

g2 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[7]][,coi] + 1),
  xlab = "Ref (log)", ylab = "Obs (log)", method = "spearman"
)

g3 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[3]]$estimate[,coi],
  xlab = "Ref", ylab = "SAVER", method = "spearman"
)

g4 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[3]]$estimate[,coi] + 1),
  xlab = "Ref (log)", ylab = "SAVER (log)", method = "spearman"
)

g5 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[4]][,coi],
  xlab = "Ref", ylab = "ALRA", method = "spearman"
)

g6 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = dat[[datoi]][[4]][,coi],
  xlab = "Ref (log)", ylab = "ALRA", method = "spearman"
)

g7 <- plot.cor(
  x = dat[[datoi]][[6]][,coi], y = dat[[datoi]][[5]][,coi],
  xlab = "Ref", ylab = "DCA", method = "spearman"
)

g8 <- plot.cor(
  x = log(dat[[datoi]][[6]][,coi] + 1), y = log(dat[[datoi]][[5]][,coi]),
  xlab = "Ref (log)", ylab = "DCA (log)", method = "spearman"
)

gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, ncol = 4)


# gene-gene correlation scatter plot
goi <- 1

g1 <- plot.cor(
  x = dat[[datoi]][[6]][goi,], y = dat[[datoi]][[7]][goi,],
  xlab = "Ref", ylab = "Obs"
)

g2 <- plot.cor(
  x = log(dat[[datoi]][[6]][goi,] + 1), y = log(dat[[datoi]][[7]][goi,] + 1),
  xlab = "Ref (log)", ylab = "Obs (log)"
)

g3 <- plot.cor(
  x = dat[[datoi]][[6]][goi,], y = dat[[datoi]][[3]]$estimate[goi,],
  xlab = "Ref", ylab = "SAVER"
)

g4 <- plot.cor(
  x = log(dat[[datoi]][[6]][goi,] + 1), y = log(dat[[datoi]][[3]]$estimate[goi,] + 1),
  xlab = "Ref (log)", ylab = "SAVER (log)"
)

g5 <- plot.cor(
  x = dat[[datoi]][[6]][goi,], y = dat[[datoi]][[4]][goi,],
  xlab = "Ref", ylab = "ALRA"
)

g6 <- plot.cor(
  x = log(dat[[datoi]][[6]][goi,] + 1), y = dat[[datoi]][[4]][goi,],
  xlab = "Ref (log)", ylab = "ALRA"
)

g7 <- plot.cor(
  x = dat[[datoi]][[6]][goi,], y = dat[[datoi]][[5]][goi,],
  xlab = "Ref", ylab = "DCA"
)

g8 <- plot.cor(
  x = log(dat[[datoi]][[6]][goi,] + 1), y = log(dat[[datoi]][[5]][goi,]),
  xlab = "Ref (log)", ylab = "DCA (log)"
)

gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, ncol = 4)


# boxplot for correlation results
cor.dat.mat <- list("gene" = vector("list",4), "cell" = vector("list",4))
names(cor.dat.mat[[1]]) <- names(dat)
names(cor.dat.mat[[2]]) <- names(dat)

for(i in 1:4){
  cor.dat.mat[[1]][[i]] <- matrix(unlist(cor.dat[[1]][[i]]), nrow = 4, byrow = T)
  rownames(cor.dat.mat[[1]][[i]]) <- c("Obs","SAVER","ALRA","DCA")
  cor.dat.mat[[1]][[i]] <- melt(cor.dat.mat[[1]][[i]])
  cor.dat.mat[[1]][[i]]$Var2 <- names(dat)[i]
  
  cor.dat.mat[[2]][[i]] <- matrix(unlist(cor.dat[[2]][[i]]), nrow = 4, byrow = T)
  rownames(cor.dat.mat[[2]][[i]]) <- c("Obs","SAVER","ALRA","DCA")
  cor.dat.mat[[2]][[i]] <- melt(cor.dat.mat[[2]][[i]])
  cor.dat.mat[[2]][[i]]$Var2 <- names(dat)[i]
}

ggplot() + theme_cowplot() + 
  geom_boxplot(
    data = rbind(cor.dat.mat[[1]][[1]],cor.dat.mat[[1]][[2]],cor.dat.mat[[1]][[3]],cor.dat.mat[[1]][[4]]),
    aes(x = Var2, y = value, col = Var1)
  )

ggplot() + theme_cowplot() + 
  geom_boxplot(
    data = rbind(cor.dat.mat[[2]][[1]],cor.dat.mat[[2]][[2]],cor.dat.mat[[2]][[3]],cor.dat.mat[[2]][[4]]),
    aes(x = Var2, y = value, col = Var1)
  )




###############################################################################
## Correlation with random reference plots

# X is ref, cell in cols
get.rand.cor.cell <- function(X, Y, num = 100, method = "spearman"){
  n.X <- ncol(X)
  rand.X <- sample(n.X, size = num)
  
  X.cor <- cor(X[,rand.X], Y, method = method)
  
  # remove entry with self
  for(ii in 1:num){
    X.cor[ii,rand.X[ii]] <- NA
  }
  
  return(X.cor)
  
  res <- colMeans(X.cor, na.rm = T)
  
  
  # for(ii in 1:n.X){
  #   if(ii %in% rand.X){
  #     rand.X.ii <- rand.X[-which(rand.X == ii)]
  #   } else {
  #     rand.X.ii <- rand.X
  #   }
  #   res[ii] <- mean(sapply(rand.X.ii, function(xx){
  #     cor(X[,xx], Y[,ii], method = method)
  #   }))
  # }
  
  return(res)
}

cor.rand.dat <- vector("list", 2)
names(cor.rand.dat) <- c("gene", "cell")

for (i in 1:2) {
  cor.rand.dat[[i]] <- vector("list", 4)
  names(cor.rand.dat[[i]]) <- names(dat)
  for (j in 1:4) {
    cor.rand.dat[[i]][[j]] <- vector("list", 4)
    names(cor.rand.dat[[i]][[j]]) <- c("Obs","SAVER","ALRA","DCA")
  }
}

n.random <- floor(n.cells / 5)
for(i in 1:4){
  for(j in 1:4){
    ind <- c(7,3,4,5)
    if(j != 2){
      cor.rand.dat[[2]][[i]][[j]] <- get.rand.cor.cell(
        X = dat[[i]][[6]], Y = dat[[i]][[ind[j]]], num = n.random[i]
      )
    } else {
      cor.rand.dat[[2]][[i]][[j]] <- get.rand.cor.cell(
        X = dat[[i]][[6]], Y = dat[[i]][[ind[j]]]$estimate, num = n.random[i]
      )
    }
  }
}

cor.norm.dat <- cor.dat
for(i in 1:4){
  for(j in 1:4){
    cor.norm.dat[[2]][[i]][[j]] <- cor.dat[[2]][[i]][[j]] / 
      apply(cor.rand.dat[[2]][[i]][[j]], 2, function(x) mean(abs(x), na.rm = T))
    # cor.norm.dat[[2]][[i]][[j]] <- cor.dat[[2]][[i]][[j]] / apply(cor.rand.dat[[2]][[i]][[j]], 2, max, na.rm = T)
  }
}


# boxplot for correlation results
cor.norm.dat.mat <- list("gene" = vector("list",4), "cell" = vector("list",4))
names(cor.norm.dat.mat[[1]]) <- names(dat)
names(cor.norm.dat.mat[[2]]) <- names(dat)

for(i in 1:4){
  cor.norm.dat.mat[[2]][[i]] <- matrix(unlist(cor.norm.dat[[2]][[i]]), nrow = 4, byrow = T)
  rownames(cor.norm.dat.mat[[2]][[i]]) <- c("Obs","SAVER","ALRA","DCA")
  cor.norm.dat.mat[[2]][[i]] <- melt(cor.norm.dat.mat[[2]][[i]])
  cor.norm.dat.mat[[2]][[i]]$Var2 <- names(dat)[i]
}

ggplot() + theme_cowplot() + 
  geom_boxplot(
    data = rbind(cor.norm.dat.mat[[2]][[1]],cor.norm.dat.mat[[2]][[2]],cor.norm.dat.mat[[2]][[3]],cor.norm.dat.mat[[2]][[4]]),
    aes(x = Var2, y = value, col = Var1)
  )


