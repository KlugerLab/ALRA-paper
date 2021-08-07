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
library(fastRPCA)
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
set.seed(3)
# subsample from PBMC library size
simCellSize <- sample(PBMC_lib_size, size = nCell, replace = F)

# get prob


bulkProb <- sweep(bulkCount, MARGIN = 2, STATS = colSums(bulkCount), FUN = "/")
cellType <- sample(bulkType, size = nCell, replace = T)

# simulate cell type
true.mat <- bulkCount[,cellType]
true.p <- sweep(true.mat,2, colSums(true.mat), '/')
expected.mat <- sweep(true.p,2, simCellSize, '*')

keep.genes <- which ( rowSums(expected.mat) > 25)
keep.cells <- which ( colSums(expected.mat) > 1000)
length(keep.genes)
true.p.sub <- true.p[keep.genes,keep.cells]

set.seed(3)
gois <- sample(names(which(rowSums(true.p.sub ==0 ) >0  )), 100,replace=F)

cells_gois_list <- list()

print("begin its")
its <- 500
for (it in 1:its)
{
    print(it)
    cells_gois <- matrix(0, length(gois), ncol(true.p.sub))
    # get count from multinomial distribution
    # simulation data matrix
    simCount <- matrix(0, nrow = nrow(true.p.sub), ncol = ncol(true.p.sub))
    for(i in 1:ncol(true.p.sub)){
        simCount[,i] <- rmultinom(n = 1, size = simCellSize[i], prob = true.p.sub[,i])
    }
    rownames(simCount)             <- rownames(true.p.sub)
    colnames(simCount)             <- colnames(true.p.sub)
    #simCount <- simCount[rowSums(simCount) > 0,]
    simNorm <- log(normbycol(simCount, 10000) + 1)
    k <- 10
    svdout <- fastPCA(simNorm,k,its=10, l=k+10)
    lra <- svdout$U %*% svdout$S %*% t( svdout$V )
    rownames(lra) <- rownames(simNorm)
    cells_gois_list[[it]] <- lra[gois,]
}


print("begin genes")

gi <- 1
cormat.means <- rep(-1,length(gois))
gplots <- list()
for (gi in 1:length(gois)) {
    print(gi)
    gene.across.its <- matrix(-1,its, ncol(true.p.sub))
    for (it in 1:its) {
       gene.across.its[it,] <- cells_gois_list[[it]][gois[gi],] 
    }
    biozeros <- which(true.p.sub[gois[gi],]== 0)
    cormat <- cor(gene.across.its[,biozeros],gene.across.its[,biozeros])
    toplot <- data.frame(  x=c(cormat- diag(rep(1,nrow(cormat)))))
    if (gi <4) {
        g1 <- ggplot(toplot) + geom_histogram(aes(x=x)) + xlim(-1,1) + ggtitle(gois[gi])+theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
        gplots[[gi]] <- g1
    }
    cormat.means[gi] <- mean(toplot$x)
}
g1

g<- gridExtra::grid.arrange(gplots[[1]], gplots[[2]], gplots[[3]], nrow=1)
ggsave('figs/simulation_cors.pdf',g)
plot(g)


print(mean(cormat.means))
print(sd(cormat.means))

