library(reldist)
library(ggplot2)
library(SAVER)
library(gridExtra)
source('/Users/george/Research_Local/ALRA/alra.R')
source('convenience.R')

if (!file.exists('data/melanoma_dropseq.rds')) {
	dropseq <- t(read.delim('data/GSE99330_dropseq_counts.txt',sep="\t"))
	saveRDS(dropseq, "data/melanoma_dropseq.rds")
}else{
	print("Reading existing dropseq file")
	dropseq <- readRDS("data/melanoma_dropseq.rds")
}
fish <- data.frame(read.delim('data/fishSubset.txt',sep=" "))


set.seed(3)
n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- colnames(dropseq)
cell.names <- rownames(dropseq)
genes <- which(gene.names %in% colnames(fish))

################################
# Run SAVER and SLRA
###############################
# We do not exclude any cells. Run SAVER on the genes of interest. 
if (!file.exists('data/melanoma_dropseq_saver.rds')) {
	system.time(saver_out <- saver(t(dropseq),pred.genes = which( colnames(dropseq) %in% ginigenes),pred.genes.only = TRUE, do.fast=F))
	saveRDS(saver_out, "data/melanoma_dropseq_saver.rds")
}else{
	print("Reading existing saver output")
	saver_out <- readRDS( "data/melanoma_dropseq_saver.rds")
}
set.seed(3)
lambda.samp <- sample.saver(saver_out)


# Compute SLRA
dropseq_norm <- normalize_data((dropseq))
k_choice <- choose_k(dropseq_norm)
k <- k_choice$k
print(sprintf("Chose k=%d",k))
system.time(result.completed <- alra(dropseq_norm,k))
slra <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))


fish.norm <- fish
dropseq.norm <- t(dropseq[,genes])
lambda.samp.norm <- lambda.samp
slra.norm <- slra[genes,]

fish.gini <- apply(fish.norm[, gene.names[genes]], 2, function(x)  gini(x[complete.cases(x)]))
dropseq.gini <- apply(dropseq.norm, 1, gini)
saver.gini <- apply(lambda.samp.norm, 1, gini)
slra.gini <- apply(slra.norm, 1, gini)
dropseq.cor <- cor(fish.gini,dropseq.gini)
saver.cor <- cor(fish.gini,saver.gini)
slra.cor <- cor(fish.gini,slra.gini)


library(cowplot)
toplot <- data.frame(FISH=fish.gini,Dropseq=dropseq.gini, SAVER=saver.gini,SLRA=slra.gini)
get.annotation <- function(base,x) {
	base+ geom_point(size=2)  + annotate(x=0.87,y=0.09, label=sprintf("r = %.2f", x),geom="text")+ xlim(0,1) + ylim(0,1)  + geom_abline()  + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
}
g1 <- ggplot(toplot, aes(x=FISH, y=Dropseq))  
g1 <- get.annotation(g1,dropseq.cor)  
g2 <- ggplot(toplot, aes(x=FISH, y=SAVER))
g2 <- get.annotation(g2,saver.cor)
g3 <- ggplot(toplot, aes(x=FISH, y=SLRA)) 
g3 <- get.annotation(g3,slra.cor)
g <-arrangeGrob(g1,g3,g2,nrow=1)
ggsave('figures/fish_comparison_saver_all.pdf',width = 7.2, height=2.,units = "in",g)
