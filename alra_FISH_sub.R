# The down-sample analysis was inspired by Huang et al.                                   
# (2018). In order to compare with their results, we closely followed their               
# analysis, the code for which is available at                                            
# https://github.com/mohuangx/SAVER-paper, and data available at                          
# https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0 


library(reldist)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(SAVER)
library(Rmagic)
library(scImpute)
source('../ALRA/alra.R')


################################
# Load and filter
###############################
if (!file.exists('data/melanoma_dropseq_filtered.rds')) {
  melanoma.raw <- as.matrix(read.table("data/GSE99330_dropseqUPM.txt"))
  # convert upm to counts
  melanoma <- sweep(melanoma.raw, 2, apply(melanoma.raw, 2, function(x) min(x[x!= 0])), "/")
  melanoma <- round(melanoma)
  # filter data
  dropseq <- melanoma[which(rowMeans(melanoma) > 0.01),
                      which(colSums(melanoma) >= 500 &
                              colSums(melanoma) <= 20000)]
  rm(melanoma.raw)
  gc()
  write.csv(dropseq,'data/GSE99330_counts_filtered.csv' )
  saveRDS(dropseq,  "data/melanoma_dropseq_filtered.rds")
}else{
  print("Loading existing melanoma data")
  dropseq <- readRDS("data/melanoma_dropseq_filtered.rds")
}


# FISH data
fish <- read.table("data/fishSubset.txt", header = TRUE, row.names = 1)


# basics
n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- rownames(dropseq)
cell.names <- colnames(dropseq)
genes <- which(gene.names %in% colnames(fish))



################################
# Run Imputation
###############################
if (!file.exists('data/melanoma_dropseq_saver.rds')) {
  system.time(saver<- saver((dropseq),pred.genes = which( rownames(dropseq) %in% colnames(fish)),pred.genes.only = TRUE, do.fast=F))
  saveRDS(saver, "data/melanoma_dropseq_saver.rds")
}else{
  print("Loading existing saver output")
  saver<- readRDS( "data/melanoma_dropseq_saver.rds")
}
set.seed(3)
saver.mat <- sample.saver(saver)

# Compute ALRA
dropseq_t <- t(dropseq)
dropseq_norm <- normalize_data(dropseq_t)
k_choice <- choose_k(dropseq_norm)

(k <- k_choice$k)
df <- data.frame(x=2:100,y=diff(k_choice$d))[4:99,]
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_vline(xintercept=k+1)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) 
ggsave('figs/spectrum_fish.pdf', width=4,height=2, g)

print(sprintf("Chose k=%d",k))
system.time(result.completed <- alra(dropseq_norm,k))
alra.mat <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))

# Compute DCA
fn <- "data/dca_fish/mean_norm.tsv"
if (!file.exists(fn)) {
    system.time(system(sprintf('dca data/GSE99330_counts_filtered.csv data/dca_fish')))
}
x.dca <- read.table(fn, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
dca.mat <- log(as.matrix(x.dca[gene.names,cell.names]) + 1)

# MAGIC
fn ='data/fish_magic.RDS';
if ( !file.exists(fn)) {
        print(sprintf("Generating %s\n", fn));
        system.time(    magic_result <- magic(dropseq_norm))
        saveRDS(magic_result,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        magic_result <- readRDS(fn)
}
magic.mat <- t(magic_result$result)
rownames(magic.mat) <- colnames(dropseq_norm)


# scImpute
set.seed(3)
fno =  "data/fish_scimpute_out/scimpute_count.csv"
if ( !file.exists(fno)) {
    print(sprintf("Generating %s\n", fno));
    fni <- "data/GSE99330_counts_filtered.csv"
    scimpute_result <- scimpute( count_path =fni, out_dir = "data/fish_scimpute_out/",
    Kcluster=5)
}

scimpute.mat <- as.matrix(read.table( fno,
    sep = ",",
    header = T,
    row.names = 1,
    stringsAsFactors = F))





################################
# Filter and Normalize
###############################

# normalize function
normalize.gapdh <- function(x) {
  x.filt <- x[, x["GAPDH", ] < quantile(x["GAPDH", ], 0.9) & x["GAPDH", ] > quantile(x["GAPDH", ], 0.1)]
  x.norm <- sweep(x.filt, 2, x.filt["GAPDH", ]/mean(x.filt["GAPDH", ]), "/")
}


# filter
fish.filt <- fish[fish[, "GAPDH"] < quantile(fish[, "GAPDH"], 0.9) & 
                    fish[, "GAPDH"] > quantile(fish[, "GAPDH"], 0.1), ]
dropseq.filt <- dropseq[genes, dropseq["GAPDH", ] < quantile(dropseq["GAPDH", ], 0.9) &
                          dropseq["GAPDH", ] > quantile(dropseq["GAPDH", ], 0.1)]
saver.mat.filt <- saver.mat[, saver.mat["GAPDH", ] < quantile(saver.mat["GAPDH", ], 0.9) &
                                  saver.mat["GAPDH", ] > quantile(saver.mat["GAPDH", ], 0.1)]
alra.filt <- alra.mat[, alra.mat["GAPDH", ] < quantile(alra.mat["GAPDH", ], 0.9) & alra.mat["GAPDH", ] > quantile(alra.mat["GAPDH", ], 0.1)]
dca.mat.filt <- dca.mat[genes, dca.mat["GAPDH", ] < quantile(dca.mat["GAPDH", ], 0.9) &
                       dca.mat["GAPDH", ] > quantile(dca.mat["GAPDH", ], 0.1)]
magic.filt <- magic.mat[, magic.mat["GAPDH", ] < quantile(magic.mat["GAPDH", ], 0.9) &
                      magic.mat["GAPDH", ] > quantile(magic.mat["GAPDH", ], 0.1)]
scimpute.filt <- scimpute.mat[, scimpute.mat["GAPDH", ] < quantile(scimpute.mat["GAPDH", ], 0.9) & scimpute.mat["GAPDH", ] > quantile(scimpute.mat["GAPDH", ], 0.1)]


# normalize
fish.norm <- sweep(fish.filt, 1, fish.filt[, "GAPDH"]/mean(fish.filt[, "GAPDH"]), "/")
dropseq.norm <- normalize.gapdh(dropseq[genes, ])
alra.norm <- normalize.gapdh(alra.mat[genes,])
saver.mat.norm <- sweep(saver.mat.filt, 2, saver.mat.filt["GAPDH", ]/mean(saver.mat.filt["GAPDH", ]), "/")
dca.norm <- normalize.gapdh(dca.mat[genes,])
magic.norm <- normalize.gapdh(magic.mat[genes,])
scimpute.norm <- normalize.gapdh(scimpute.mat[genes,])




## Gini

# calculate gini
fish.gini <- apply(fish.norm[, gene.names[genes][-6]], 2, function(x)  gini(x[complete.cases(x)]))
dropseq.gini <- apply(dropseq.norm[-6,], 1, gini)
saver.gini <- apply(saver.mat.norm[-6,], 1, gini)
alra.gini <- apply(alra.norm[-6,], 1, gini)
dca.gini <- apply(dca.norm[-6,], 1, gini)
magic.gini <- apply(magic.norm[-6,], 1, gini)
scimpute.gini <- apply(scimpute.norm[-6,], 1, gini)


# cor
dropseq.cor <- cor(fish.gini,dropseq.gini)
saver.cor <- cor(fish.gini,saver.gini)
alra.cor <- cor(fish.gini,alra.gini)
dca.cor <- cor(fish.gini,dca.gini)
magic.cor <- cor(fish.gini,magic.gini)
scimpute.cor <- cor(fish.gini,scimpute.gini)


# plot
toplot <- data.frame(
  FISH=fish.gini, Dropseq=dropseq.gini, SAVER=saver.gini,ALRA=alra.gini, DCA=dca.gini, MAGIC=magic.gini, scImpute=scimpute.gini
)

get.annotation <- function(base,x) {
  base + geom_point(size=2) + theme_cowplot() + 
    annotate(x=0.87,y=0.09, label=sprintf("r = %.2f", x),geom="text")+ xlim(0,1) + ylim(0,1) + 
    geom_abline()  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
}

g1 <- ggplot(toplot, aes(x=FISH, y=Dropseq)) + ggtitle('Original')
g1 <- get.annotation(g1,dropseq.cor)
g2 <- ggplot(toplot, aes(x=FISH, y=ALRA))+ ggtitle('ALRA')
g2 <- get.annotation(g2,alra.cor)
g3 <- ggplot(toplot, aes(x=FISH, y=SAVER)) + ggtitle('SAVER')
g3 <- get.annotation(g3,saver.cor) 
g4 <- ggplot(toplot, aes(x=FISH, y=DCA))+ ggtitle('DCA')
g4 <- get.annotation(g4,dca.cor) 
g5 <- ggplot(toplot, aes(x=FISH, y=MAGIC))+ ggtitle('MAGIC')
g5 <- get.annotation(g5,magic.cor) 
g6 <- ggplot(toplot, aes(x=FISH, y=scImpute))+ ggtitle('scImpute')
g6 <- get.annotation(g6, scimpute.cor)
g <- gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6, nrow=2)

ggsave( filename = "figs/fish_comparison_sub.pdf",
       plot = g,
       width = 9, height = 5, units = "in")


