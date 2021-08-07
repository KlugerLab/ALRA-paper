library(reldist)
library(ggplot2)
library(SAVER)
library(gridExtra)
source('../ALRA/alra.R')


## Load data

#DAta can be obtained from: https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0
   if (!file.exists('data/melanoma_dropseq.rds')) {
    df.upm <- as.matrix(read.table("data/GSE99330_dropseqUPM.txt", row.names = 1, header = TRUE))
    dropseq <- sweep(df.upm, 2, apply(df.upm, 2, function(x) min(x[x!= 0])), "/")
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
   sum(!is.wholenumber(dropseq))
   dropseq <- round(dropseq)  
   write.csv(dropseq,'data/GSE99330_counts.csv' )
   #dropseq <- t(dropseq)
  saveRDS(dropseq, "data/melanoma_dropseq.rds")
}else{
  print("Reading existing dropseq file")
  dropseq <- readRDS("data/melanoma_dropseq.rds")
}
# FISH data
fish <- read.table("data/fishSubset.txt", header = TRUE, row.names = 1)


# basics
set.seed(3)
n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- rownames(dropseq)
cell.names <- colnames(dropseq)
genes <- which(gene.names %in% colnames(fish))


which( colnames(dropseq) %in% colnames(fish))
colnames(fish)%in% rownames(dropseq)
################################
# Run Imputation
###############################
sum(rowSums(dropseq) ==0)
which( rownames(dropseq) %in% colnames(fish))

# We do not exclude any cells. Run SAVER on the genes of interest. 
if (!file.exists('data/melanoma_dropseq_saver.rds')) {
  system.time(saver_out <- saver((dropseq),pred.genes =which(rownames(dropseq)%in% colnames(fish)),pred.genes.only = TRUE, do.fast=F))
  saveRDS(saver_out, "data/melanoma_dropseq_saver.rds")
}else{
  print("Reading existing saver output")
  saver_out <- readRDS( "data/melanoma_dropseq_saver.rds")
}
set.seed(3)
lambda.samp <- sample.saver(saver_out)


# Compute ALRA
dropseq_norm <- normalize_data(t(dropseq))
k_choice <- choose_k(dropseq_norm)
k <- k_choice$k
print(sprintf("Chose k=%d",k))
system.time(result.completed <- alra(dropseq_norm,k))
alra.mat <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))


# DCA results
# read in dca result
fn <- "data/dca_fish/mean_norm.tsv"
if (!file.exists(fn)) {
    system.time(system(sprintf('dca data/GSE99330_counts.csv data/dca_fish')))
}
x.dca <- read.table(fn, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
dca_norm <- log(as.matrix(x.dca[gene.names,cell.names]) + 1)

## Filter & Normalize

# normalize
dropseq.norm <- dropseq[genes,]
fish.norm <- fish
lambda.samp.norm <- lambda.samp
alra.mat.norm <- alra.mat[genes,]
dca.norm <- dca_norm[genes,]



## Gini

# calculate gini
fish.gini <- apply(fish.norm[, gene.names[genes]], 2, function(x)  gini(x[complete.cases(x)]))
dropseq.gini <- apply(dropseq.norm, 1, gini)
saver.gini <- apply(lambda.samp.norm, 1, gini)
alra.mat.gini <- apply(alra.mat.norm, 1, gini)
dca.gini <- apply(dca.norm, 1, gini)


# cor
dropseq.cor <- cor(fish.gini,dropseq.gini)
saver.cor <- cor(fish.gini,saver.gini)
alra.mat.cor <- cor(fish.gini,alra.mat.gini)
dca.cor <- cor(fish.gini,dca.gini)


# plot
toplot <- data.frame(FISH=fish.gini, Dropseq=dropseq.gini, SAVER=saver.gini, SLRA=alra.mat.gini, DCA=dca.gini)

get.annotation <- function(base,x) {
  base + geom_point(size=2) + 
    annotate(x=0.87,y=0.09, label=sprintf("r = %.2f", x),geom="text")+ xlim(0,1) + ylim(0,1) + 
    geom_abline()  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
}

g1 <- ggplot(toplot, aes(x=FISH, y=Dropseq))
g1 <- get.annotation(g1,dropseq.cor)

g2 <- ggplot(toplot, aes(x=FISH, y=SAVER))
g2 <- get.annotation(g2,saver.cor)

g3 <- ggplot(toplot, aes(x=FISH, y=SLRA)) 
g3 <- get.annotation(g3,alra.mat.cor)

g4 <- ggplot(toplot, aes(x=FISH, y=DCA))
g4 <- get.annotation(g4,dca.cor)

ggsave("./fish_logdca.pdf", width = 2.4, height=2., units = "in", plot = g4)


