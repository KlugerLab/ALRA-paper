# The analysis of Hrvatin et al. (2018) data was inspired by Huang et al.
# (2018). In order to compare with their results, we closely followed their
# analysis, the code for which is available at
# https://github.com/mohuangx/SAVER-paper, and data available at
# https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0
source('convenience.R')
library(gridExtra)
library(ggplot2)
library(Seurat)
library(ggplot2)
library(cowplot)
source('../ALRA/alra.R')

if (!file.exists('data/hrvatin_full.rds')) {
	x <- read.csv("SAVER-data/GSE102827_merged_all_raw.csv", header = TRUE, row.names = 1, check.names = FALSE)
	x <- as.matrix(x)
	# Filter out genes with expression less than 0.00003
	x1 <- x[rowMeans(x) >= 0.00003, ]
	# look at non-zero expression
	# Filter out genes with non-zero expression in less than 4 cells
	x2 <- x1[rowSums(x1 != 0) >= 4, ]
	dim(x2)
	#set.seed(011118)
	set.seed(3)
	samp.cells <- sample(1:ncol(x2), 20000)
	samp.cells <- 1:ncol(x2)
	x.sub <- x2[, samp.cells]
	saveRDS(x.sub, "data/hrvatin_full.rds")
	x <- x.sub
}else{
	print("Reading")
	x <- readRDS('data/hrvatin_full.rds')
}


# Compute ALRA
if (!file.exists('data/hrvatin_full_alra.rds')) {
	print("Computing ALRA")
	x_t <- t(x)
	totalUMIPerCell <- rowSums(x_t);
	x_norm <- sweep(x_t, 1, totalUMIPerCell, '/');
	x_norm <- x_norm * 10E3
	x_norm <- log(x_norm +1);

	set.seed(3)
	system.time(k_choice <- choose_k(x_norm))
	saveRDS(k_choice,file='data/k_choice.rds')


	print(sprintf("k=%d was chosen", k_choice$k))
	system.time(result.completed <- alra(x_norm,k_choice$k)) 
	x.alra <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
	sum(x.alra==0)/(nrow(x.alra)*ncol(x.alra))
	colnames(x.alra) <- colnames(x)
	saveRDS(x.alra,file='data/hrvatin_full_alra.rds')
}else{
	print("Reading ALRA")
	x.alra <- readRDS('data/hrvatin_full_alra.rds')
	k_choice <- readRDS('data/k_choice.rds')
}

df <- data.frame(x=2:100,y=diff(k_choice$d))[10:99,]
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=2)  + geom_line()+ geom_vline(xintercept=k_choice$k +1)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,90,10))
ggsave('figs/spectrum_hrvatin.pdf', width=4,height=2, g)


ct.dat <- read.csv("SAVER-data/hrvatin_celltypes.csv", header = TRUE,
	    row.names = 1, stringsAsFactors = FALSE)
ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
ct.dat2$celltype[!is.na(ct.dat2$subtype)] <- 
  ct.dat2$subtype[!is.na(ct.dat2$subtype)]
table(ct.dat2$celltype)
ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]
ident <- ct.dat3[colnames(x), 4]
x.sub <- x[, !is.na(ident)]
x.alra.sub <- x.alra[, !is.na(ident)]
ident2 <- ct.dat3[colnames(x.sub), 4]



obs <- CreateSeuratObject(x.sub, project = "obs")
alra <- CreateSeuratObject(x.alra.sub, project = "alra")

obs <- NormalizeData(obs)

obs <- FindVariableGenes(obs, do.plot = FALSE)
alra <- FindVariableGenes(alra, do.plot = FALSE)

obs <- ScaleData(obs)
alra <- ScaleData(alra)

obs <- RunPCA(obs, pc.genes = obs@var.genes, do.print = FALSE, pcs.compute = 70)
alra <- RunPCA(alra, pc.genes = alra@var.genes, do.print = F, pcs.compute = 70)


ident3 <- factor(ident2, levels = levels(ident2)[c(3, 4, 5, 6, 7, 8, 9, 1, 2, 10, 12, 13, 11)])
dim.obs <- k_choice$k
dim.alra <- k_choice$k
obs <- RunTSNE(obs, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", fast_tsne_path="/data/george/Research_Local/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12)
alra <- RunTSNE(alra, dims.use = 1:dim.alra, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", fast_tsne_path="/data/george/Research_Local/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12)
obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
alra.tsne <- data.frame(alra@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)


print("Done with t-SNE, beginning random forests")
library(ranger)
df <- data.frame(t(obs@scale.data[obs@var.genes,]), Y=as.factor(ident2))
rfr <- ranger(Y ~ ., df)
errs <- df$Y != rfr$predictions
df <- data.frame(t(x.alra.sub[alra@var.genes,]), Y=as.factor(ident2))
rfr.alra <- ranger(Y ~ ., df)
errs.alra <- df$Y != rfr.alra$predictions
rfr
rfr.alra





library(ggplot2)
library(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my.cols <- gg_color_hue(length(unique(ident2)))
p3 <- ggplot(obs.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
	    size = 0.2) + guides(colour = FALSE) + 
	theme(axis.title = element_blank(), axis.text = element_blank(), 
	      axis.ticks = element_blank(), 
	      plot.title = element_text(size = 16, face = "plain")) + 
		labs(title = "Observed") +annotate(x=max(obs.tsne$tSNE_1)-10,
		y=min(obs.tsne$tSNE_2), 
		label=sprintf("OOB Error = %.1f%%", rfr$prediction.error*100),
		size=4,geom="text")
p2 <- ggplot(alra.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
	size = 0.2)  + theme(axis.title = element_blank(), axis.text = element_blank(), 
	axis.ticks = element_blank(), 
	plot.title = element_text(size = 16, face = "plain"), 
	legend.position="bottom", legend.title=element_blank()) + labs(title = "ALRA")+
	scale_colour_manual(name = "Cell type", values = my.cols)  + 
	guides(colour=guide_legend(override.aes=list(size=8), nrow=4))+
	annotate(x=max(alra.tsne$tSNE_1)-10,y=min(alra.tsne$tSNE_2), 
	 label=sprintf("OOB Error = %.1f%%", rfr.alra$prediction.error*100),
	 size=4,geom="text")
g1<-plot_grid(p3,p2+theme(legend.position="none"),nrow=1)
ggsave(g1,filename="figs/Hrvatin-tsne.pdf",width=9,height=4.5)

gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}
g2 <- gglegend(p2)
ggsave(g2,filename="figs/Hrvatin-tsne_legend.pdf",width=12,height=2)
