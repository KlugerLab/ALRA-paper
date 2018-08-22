##################################
##  Using same distribution
#################################
library(SAVER)
source('convenience.R')
library(reldist)
library(gridExtra)
library(ggplot2)
source('../ALRA/alra.R')


################################
# Load and filter
###############################
if (!file.exists('data/melanoma_dropseq_filtered.rds')) {
	melanoma.raw <- as.matrix(read.table("data/GSE99330_dropseqUPM.txt.gz"))
	# convert upm to counts
	melanoma <- sweep(melanoma.raw, 2, apply(melanoma.raw, 2,
						 function(x) min(x[x!= 0])), "/")
	melanoma <- round(melanoma)
	# filter data
	dropseq <- melanoma[which(rowMeans(melanoma) > 0.01),
				  which(colSums(melanoma) >= 500 &
					  colSums(melanoma) <= 20000)]
	rm(melanoma.raw)
	gc()
	saveRDS(dropseq,  "data/melanoma_dropseq_filtered.rds")
}else{
	print("Loading existing melanoma data")
	dropseq <- readRDS("data/melanoma_dropseq_filtered.rds")
}


fish <- read.table("data/fishSubset.txt", header = TRUE, row.names = 1)

n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- rownames(dropseq)
cell.names <- colnames(dropseq)
genes <- which(gene.names %in% colnames(fish))



################################
# Run SAVER and SLRA
###############################
if (!file.exists('data/melanoma_dropseq_filtered_saver.rds')) {
	system.time(saver<- saver((dropseq),pred.genes = which( rownames(dropseq) %in% colnames(fish)),pred.genes.only = TRUE, do.fast=F))
	saveRDS(saver, "data/melanoma_dropseq_filtered_saver.rds")
}else{
	print("Loading existing saver output")
	saver<- readRDS( "data/melanoma_dropseq_filtered_saver.rds")
}
set.seed(3)
lambda.samp <- sample.saver(saver)

# Compute SLRA
dropseq_t <- t(dropseq)
dropseq_norm <- normalize_data(dropseq_t)
k_choice <- choose_k(dropseq_norm)

k <- k_choice$k
df <- data.frame(x=2:100,y=diff(k_choice$d))[4:99,]
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_vline(xintercept=k+1)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) 
ggsave('figures/spectrum_fish.pdf', width=4,height=2, g)

print("Chose k=%d",k)
system.time(result.completed <- alra(dropseq_norm,k))
slra <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))

################################
# Filter and Normalize
###############################

normalize.gapdh <- function(x) {
  x.filt <- x[, x["GAPDH", ] < quantile(x["GAPDH", ], 0.9) & x["GAPDH", ] > quantile(x["GAPDH", ], 0.1)]
  x.norm <- sweep(x.filt, 2, x.filt["GAPDH", ]/mean(x.filt["GAPDH", ]), "/")
}


fish.filt <- fish[fish[, "GAPDH"] < quantile(fish[, "GAPDH"], 0.9) & 
                    fish[, "GAPDH"] > quantile(fish[, "GAPDH"], 0.1), ]
dropseq.filt <- dropseq[genes, dropseq["GAPDH", ] < quantile(dropseq["GAPDH", ], 0.9) &
                          dropseq["GAPDH", ] > quantile(dropseq["GAPDH", ], 0.1)]
saver.filt <- saver$estimate[, saver$estimate["GAPDH", ] < 
                               quantile(saver$estimate["GAPDH", ], 0.9) & saver$estimate["GAPDH", ] >
                               quantile(saver$estimate["GAPDH", ], 0.1)]
lambda.samp.filt <- lambda.samp[, lambda.samp["GAPDH", ] < quantile(lambda.samp["GAPDH", ], 0.9) &
                    lambda.samp["GAPDH", ] > quantile(lambda.samp["GAPDH", ], 0.1)]
slra.filt <- slra[, slra["GAPDH", ] < quantile(slra["GAPDH", ], 0.9) & slra["GAPDH", ] >
                               quantile(slra["GAPDH", ], 0.1)]

fish.norm <-   sweep(fish.filt, 1, fish.filt[, "GAPDH"]/mean(fish.filt[, "GAPDH"]), "/")
dropseq.norm <- normalize.gapdh(dropseq[genes, ])
slra.norm <- normalize.gapdh(slra[genes,])
lambda.samp.norm <- sweep(lambda.samp.filt, 2, lambda.samp.filt["GAPDH", ]/mean(lambda.samp.filt["GAPDH", ]), "/")


##################
# Call GINI
###############
fish.gini <- apply(fish.norm[, gene.names[genes][-6]], 2, function(x)  gini(x[complete.cases(x)]))
dropseq.gini <- apply(dropseq.norm[-6,], 1, gini)
saver.gini <- apply(lambda.samp.norm[-6,], 1, gini)
slra.gini <- apply(slra.norm[-6,], 1, gini)
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
ggsave('figures/fish_comparision_saver_sub.pdf',width = 7.2, height=2.,units = "in",g)

