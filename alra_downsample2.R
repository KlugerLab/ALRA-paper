# The down-sample analysis was inspired by Huang et al.
# (2018). In order to compare with their results, we closely followed their
# analysis, the code for which is available at
# https://github.com/mohuangx/SAVER-paper, and data available at
# https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0


library(Seurat) # Version 2.3.4
library(cowplot)
library(ranger)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

source('convenience.R')
set.seed(3)
source("../ALRA/alra.R")


## Prepare data

# please download data with "ref" "samp.rds" "samp_saver.rds" "samp_magic.rds" "samp_scimpute.rds" patterns

dat <- vector("list", 4)
names(dat) <- c("baron", "chen", "manno", "zeisel")


# load data
ref.files <- list.files("data", pattern = "ref", full.names = TRUE)
samp.files <- list.files("data", pattern = "samp.rds", full.names = TRUE)
saver.files <- list.files("data", pattern = "samp_saver.rds", full.names = TRUE)
magic.files <- list.files("data", pattern = "samp_magic.rds", full.names = TRUE)
scimpute.files <- list.files("data", pattern = "samp_scimpute.rds", full.names = TRUE)


normalizeData <- function(x, y = x) {
  sf <- colSums(y)/mean(colSums(y))
  return(sweep(x, 2, sf, "/"))
}

# dat is a list of 4 elements, where each of them corresponds toe ach of the
# four datasets. Each element then has 9 matrices, corresponding to the
# original data, subsampled, imputations, normalized, etc.
for (i in 1:4) {
  dat[[i]] <- vector("list", 9)
  names(dat[[i]]) <- c("X", "Y", "saver", "alra", "dca", "X.norm", "Y.norm", "magic", "scimpute")
  dat[[i]][['X']] <- readRDS(ref.files[i])
  dat[[i]][['Y']] <- readRDS(samp.files[i])
  dat[[i]][['saver']] <- readRDS(saver.files[i])
  dat[[i]][['magic']] <- readRDS(magic.files[i])
  dat[[i]][['scimpute']] <- readRDS(scimpute.files[i])
  dat[[i]][['scimpute']][is.na(dat[[i]][['scimpute']])] <- 0
}
n.cells <- sapply(dat, function(x) ncol(x[[1]]))
n.genes <- sapply(dat, function(x) nrow(x[[1]]))
cell.names <- sapply(dat, function(x) colnames(x[[1]]))
gene.names <- sapply(dat, function(x) rownames(x[[1]]))
colnames(dat[[1]][['scimpute']]) <- gsub(".","-", colnames(dat[[1]][['scimpute']]),fixed=T)
colnames(dat[[3]][['scimpute']]) <- substring(colnames(dat[[3]][['scimpute']]),2)
colnames(dat[[4]][['scimpute']]) <- substring(colnames(dat[[4]][['scimpute']]),2)
for (i in 1:4) {
  dat[[i]][['X.norm']] <- normalizeData(dat[[i]][['X']])
  dat[[i]][['Y.norm']] <- normalizeData(dat[[i]][['Y']])
  dat[[i]][['scimpute']] <- normalizeData(dat[[i]][['scimpute']])
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
  dat[[i]][['alra']] <- alraResult(dat[[i]][['Y']])
}


# magic
for(i in 1:4){
  dat[[i]][['magic']] <- t(magic(normalize_data(t(dat[[i]][['Y']])))$result)
}


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
  dat[[i]][['dca']] <- dcaResult(dat[[i]][['Y']], name = names(dat)[i], do.read = T)
}




###############################################################################
## Angles

to.eval <- c( "X", "Y.norm", "saver", "alra", "dca", "magic", "scimpute")
dat.k <- vector("list",4)
names(dat.k) <- names(dat)
for(i in 1:4){
    dat.k[[i]] <- vector("list",length(dat[[i]]))
    names(dat.k[[i]]) <- names(dat[[i]])
  #for(j in 1:6){
    print(names(dat)[i])
  for(j in to.eval){
      print(j)
        if ( j == "alra" | j == "magic"  | j=="Y.norm" | j == "X") {
            data.in <- dat[[i]][[j]]
      }else if (j == "saver"){
            data.in <- dat[[i]][[j]]$estimate
        }else {
            data.in <-log(1+dat[[i]][[j]]) 
        }
        X_temp <- CreateSeuratObject(raw.data = data.in)
        X_temp@data <- data.in
        if (j == "X" | j == "saver"){
            X_temp <- NormalizeData(X_temp)
        }
        X_temp <- ScaleData(X_temp)
        X_temp <- FindVariableGenes(X_temp, do.plot = F)
        X_temp <- RunPCA(X_temp, pcs.compute = 50, do.print = F)
        original_pcs <- X_temp@dr$pca@cell.embeddings
        X_temp <- JackStraw(X_temp, num.pc = 50, prop.freq = 0.1)
        X_temp.js <- JackStrawPlot(X_temp, PCs=1:50)
        #Note that the output of JackStrawPlot for older versions of Seurat is
        #ggplot, which breaks this code. This code is tested for Seurat 2.3.4
        (chosen.k <- max(which( X_temp.js@dr$pca@jackstraw@overall.p.values[,2] < 1E-5)))
        print(X_temp.js@dr$pca@jackstraw@overall.p.values[,2]) 
        print(chosen.k)
        dat.k[[i]][[j]]  <- list( k =chosen.k, seurat = X_temp)
  }
}


dat.angles <- vector("list",4)
names(dat.angles) <- names(dat)
for(i in 1:4){
    dat.angles[[i]] <- vector("list",length(dat[[i]]))
    names(dat.angles[[i]]) <- names(dat[[i]])
    print(names(dat)[i])
  use_k <- dat.k[[i]][["X"]]$k
  original_pcs <- dat.k[[i]][["X"]]$seurat@dr$pca@cell.embeddings[,1:use_k]
  for(j in to.eval){
      imputed_pcs <- dat.k[[i]][[j]]$seurat@dr$pca@cell.embeddings[,1:use_k]
      dat.angles[[i]][j] <- pracma::subspace(as.matrix(original_pcs),as.matrix(imputed_pcs))
  }
}
print(dat.angles)
print(use_k)
dat.k[['baron']][["X"]]$k
dat.k[['chen']][["X"]]$k
dat.k[['manno']][["X"]]$k
dat.k[['zeisel']][["X"]]$k



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
  return(list(X.S@dr$tsne@cell.embeddings, X.S))
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
seurat.obj <- tsne

for(i in 1:4){
  for(j in 1:6){
    ind <- c("Y.norm","saver","alra","dca","magic","scimpute")
    if(ind[j] == "alra" | ind[j] == "magic"){
      tsne.result <- get.tsne(X = dat[[i]][[ind[j]]], n.pc = dat.k[[i]][[ind[j]]]$k)
    } else if (ind[j] == "saver") {
      tsne.result <- get.tsne(X = dat[[i]][[ind[j]]]$estimate, n.pc = dat.k[[i]][[ind[j]]]$k, do.norm = T)
    } else {
      tsne.result <- get.tsne(X = log(dat[[i]][[ind[j]]] + 1), n.pc = dat.k[[i]][[ind[j]]]$k)
    }
    tsne[[i]][[j]] <- tsne.result[[1]]
    seurat.obj[[i]][[j]] <- tsne.result[[2]]
  }
}
for(i in 1:4){
  rownames(tsne[[i]][[4]]) <- colnames(dat[[i]][[1]])
  rownames(tsne[[i]][[5]]) <- colnames(dat[[i]][[1]])
}
dat.k[[i]][ind[j]]$k



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
    #ind <- c(7,3,4,5,8,9, 1)
    ind <- c("Y.norm","saver","alra","dca","magic","scimpute" ,"X.norm")
    if(ind[j] == 'X.norm'){
      vargenes[[i]][[j]] <- get.vargenes(X = log(dat[[i]][[ind[j]]] + 1))
    } else {
      vargenes[[i]][[j]] <- seurat.obj[[i]][[j]]@var.genes
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
    # ind <- c("Y.norm","saver","alra","dca","magic","scimpute" ,"X.norm")
    # if(ind[j] == 'alra' | ind[j] == 'magic'){
    #   dat.j <- dat[[i]][[ind[j]]]
    # } else if (ind[j] == 'saver') {
    #   dat.j <- log(dat[[i]][[ind[j]]]$estimate + 1)
    # } else {
    #   dat.j <- log(dat[[i]][[ind[j]]] + 1)
    # }
    if(j == 7){
      dat.j <- log(dat[[i]][["X.norm"]] + 1)
    } else {
      dat.j <- seurat.obj[[i]][[j]]@data
    }
    if(ind[j] == 'dca' | ind[j] == 'scimpute' | ind[j] == "magic"){
      colnames(dat.j) <- colnames(dat[[i]][[1]])
    }
    oob.error[[i]][[j]] <- get.rfr(
      exprs.mat = as.matrix(dat.j[vargenes[[i]][[j]],names(info[[2]][[i]]$Ref)]), Y = info[[2]][[i]]$Ref
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


# average oob error
for(j in 1:6){
  cat(names(oob.error[[1]])[j], ": ", sep = "")
  cat(format(mean(c(oob.error[[1]][[j]]$prediction.error*100, oob.error[[2]][[j]]$prediction.error*100,
             oob.error[[3]][[j]]$prediction.error*100, oob.error[[4]][[j]]$prediction.error*100)), digits = 3))
  cat("\n")
}



## Jaccard index
library(clusteval)
library(scales)

# clustering
res <- seq(0.4, 1.4, 0.1)
ident <- vector("list", 4)
names(ident) <- names(dat)
for(i in 1:4){
  ident[[i]] <- vector("list", length = 6)
  names(ident[[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
}

for(i in 1:4){
  for(j in 1:6){
    ident[[i]][[j]] <- vector("list", 11)
    for(k in 1:11){
      seurat.obj[[i]][[j]] <- FindClusters(seurat.obj[[i]][[j]], resolution = res[k], print.output = FALSE, save.SNN = TRUE)
      ident[[i]][[j]][[k]] <- seurat.obj[[i]][[j]]@ident
    }
  }
}

# remove bad index in Baron
for(j in 1:6){
  for(k in 1:11){
    ident[[1]][[j]][[k]] <- factor(ident[[1]][[j]][[k]][match(names(info[[2]][[1]]$Ref), seurat.obj[[1]][[1]]@cell.names)])
  }
}


# calculate index
max.jaccard <- matrix(0, nrow = 4, ncol = 6)
for(i in 1:4){
  for(j in 1:6){
    max.jaccard[i, j] <- max(sapply(ident[[i]][[j]], function(x) 
      cluster_similarity(info[[2]][[i]]$Ref, x))[1:11])
  }
}
rownames(max.jaccard) <- names(dat)
colnames(max.jaccard) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
write.table(max.jaccard, "data/downsample_jaccard.txt", sep = "\t", quote = F, col.names = T, row.names = T)



## SVD
svdResult <- function(x){
  x.norm <- normalize_data(t(x))
  x.alra <- t(alra(x.norm)[[1]])
  rownames(x.alra) <- rownames(x)
  colnames(x.alra) <- colnames(x)
  return(x.alra)
}

dat.svd <- vector("list", 4)
for(i in 1:4){
  dat.svd[[i]] <- svdResult(dat[[i]][['Y']])
}

tsne.svd <- vector("list", 4)
for(i in 1:4){
  tsne.svd[[i]] <- get.tsne(X = dat.svd[[i]], n.pc = dat.k[[i]][["X"]]$k)
}

vargenes.svd <- vector("list",4)
oob.error.svd <- vector("list",4)
for(i in 1:4){
  vargenes.svd[[i]] <- get.vargenes(X = dat.svd[[i]])
  oob.error.svd[[i]] <- get.rfr(
    exprs.mat = dat.svd[[i]][vargenes.svd[[i]],names(info[[2]][[i]]$Ref)], Y = info[[2]][[i]]$Ref
  )
}

for(i in 1:4){
  gg1 <- plot.tsne(
    x = tsne.svd[[i]][names(info[[2]][[i]]$Ref),], label = info[[2]][[i]]$Ref, 
    cols.use = tsne.cols[[i]], title = "SVD", size=0.1
  ) + annotation_compass(sprintf("Error: %.1f%%", oob.error.svd[[i]]$prediction.error*100), "SE",fontsize=9 )
  ggsave(plot=gg1,filename=sprintf('figs/%s_tsne_svd.pdf', names(dat)[[i]]),width=1.45,height=1.45)
}






## In the following sections, correlation analysis is done, 

################################################################################
### Correlation with reference plots

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
 for (j in 1:6) {
   cor.dat[[i]][[j]] <- vector("list", 6)
   names(cor.dat[[i]][[j]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
 }
}

for (i in 1:4) {
 for (j in 1:6) {
   ind <- c("Y.norm","saver","alra","dca","magic","scimpute")
   if(ind[j] != "saver"){
     cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]])
     cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]])
   } else {
     cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
     cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
   }
 }
}

#
#
#
#
#
# boxplot for correlation results
cor.dat.mat <- list("gene" = vector("list",4), "cell" = vector("list",4))
names(cor.dat.mat[[1]]) <- names(dat)
names(cor.dat.mat[[2]]) <- names(dat)

for(i in 1:4){
 cor.dat.mat[[1]][[i]] <- matrix(unlist(cor.dat[[1]][[i]]), nrow = 6, byrow = T)
 rownames(cor.dat.mat[[1]][[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
 cor.dat.mat[[1]][[i]] <- melt(cor.dat.mat[[1]][[i]])
 cor.dat.mat[[1]][[i]]$Var2 <- names(dat)[i]
 cor.dat.mat[[1]][[i]]$Var1 <- factor(
   as.character(cor.dat.mat[[1]][[i]]$Var1), levels = c("Obs","ALRA","MAGIC","SAVER","scImpute","DCA")
 )

 cor.dat.mat[[2]][[i]] <- matrix(unlist(cor.dat[[2]][[i]]), nrow = 6, byrow = T)
 rownames(cor.dat.mat[[2]][[i]]) <- c("Obs","SAVER","ALRA","DCA","MAGIC","scImpute")
 cor.dat.mat[[2]][[i]] <- melt(cor.dat.mat[[2]][[i]])
 cor.dat.mat[[2]][[i]]$Var2 <- names(dat)[i]
 cor.dat.mat[[2]][[i]]$Var1 <- factor(
   as.character(cor.dat.mat[[2]][[i]]$Var1), levels = c("Obs","ALRA","MAGIC","SAVER","scImpute","DCA")
 )
}

gg1 <- ggplot() + theme_cowplot() +
 geom_boxplot(
   data = rbind(cor.dat.mat[[1]][[1]],cor.dat.mat[[1]][[2]],cor.dat.mat[[1]][[3]],cor.dat.mat[[1]][[4]]),
   aes(x = Var2, y = value, col = Var1)
 ) + ggtitle("Gene") + xlab(NULL) + ylab("Correlation with reference") + theme(legend.position = "none")

gg2 <- ggplot() + theme_cowplot() +
 geom_boxplot(
   data = rbind(cor.dat.mat[[2]][[1]],cor.dat.mat[[2]][[2]],cor.dat.mat[[2]][[3]],cor.dat.mat[[2]][[4]]),
   aes(x = Var2, y = value, col = Var1)
 ) + ggtitle("Cell") + xlab(NULL) + ylab("Correlation with reference") + theme(legend.position = "none")

ggsave(grid.arrange(gg1,gg2, ncol = 1), filename = "figs/downsample_correlation.pdf", width = 8, height = 10)
ggsave(get_legend(gg2 + theme(legend.position = "bottom", legend.spacing.x = unit(0.3,'cm'), legend.title = element_blank()) + guides(colour = guide_legend(nrow = 1))), filename = "figs/downsample_legend.pdf", width = 6, height = 1)

