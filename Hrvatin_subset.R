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


# Read in data and result from Huang et al. paper
x <- readRDS("SAVER-data/hrvatin.rds")
x.saver <- readRDS("SAVER-data/hrvatin_saver.rds")

# Compute SLRA
print("Computing SLRA")
x_t <- t(x)
totalUMIPerCell <- rowSums(x_t);
x_norm <- sweep(x_t, 1, totalUMIPerCell, '/');
x_norm <- x_norm * 10E3
x_norm <- log(x_norm +1);

set.seed(3)
system.time(k_choice <- choose_k(x_norm))
print(sprintf("k=%d was chosen", k_choice$k))
system.time(result.completed <- alra(x_norm,k_choice$k)) 
x.alra <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
sum(x.alra==0)/(nrow(x.alra)*ncol(x.alra))
colnames(x.alra) <- rownames(x_norm)


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
x.saver.sub <- x.saver[, !is.na(ident)]
ident2 <- ct.dat3[colnames(x.sub), 4]

obs <- CreateSeuratObject(x.sub, project = "obs")
alra <- CreateSeuratObject(x.alra.sub, project = "alra")
saver <- CreateSeuratObject(x.saver.sub, project = "saver")

obs <- NormalizeData(obs)
saver <- NormalizeData(saver)
#alra <- NormalizeData(alra) # SLRA is already log normalized

obs <- FindVariableGenes(obs, do.plot = FALSE)
saver <- FindVariableGenes(saver, do.plot = FALSE)
alra <- FindVariableGenes(alra, do.plot = FALSE)


library(ranger)
df <- data.frame(t(x.sub[obs@var.genes,]), Y=as.factor(ident2))
rfr <- ranger(Y ~ ., df)
df <- data.frame(t(x.alra.sub[alra@var.genes,]), Y=as.factor(ident2))
rfr.alra <- ranger(Y ~ ., df)
df <- data.frame(t(x.saver.sub[saver@var.genes,]), Y=as.factor(ident2))
rfr.saver <- ranger(Y ~ ., df)
rfr
rfr.alra
rfr.saver

obs <- ScaleData(obs)
saver <- ScaleData(saver)
alra <- ScaleData(alra)




obs <- RunPCA(obs, pc.genes = obs@var.genes, do.print = FALSE, pcs.compute = 50)
saver <- RunPCA(saver, pc.genes = saver@var.genes, do.print = FALSE, pcs.compute = 50)
alra <- RunPCA(alra, pc.genes = alra@var.genes, do.print = F, pcs.compute = 50)

dim.obs <- 35
dim.saver <- 40
dim.alra <- k_choice$k
obs <- RunTSNE(obs, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", fast_tsne_path="/data/george/Research_Local/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12)
saver <- RunTSNE(saver, dims.use = 1:dim.saver, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", fast_tsne_path="/data/george/Research_Local/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12)
alra <- RunTSNE(alra, dims.use = 1:dim.alra, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", fast_tsne_path="/data/george/Research_Local/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12)


ident3 <- factor(ident2, levels = levels(ident2)[c(3, 4, 5, 6, 7, 8, 9, 1, 2, 10, 12, 13, 11)])
obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
alra.tsne <- data.frame(alra@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)


obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
saver.tsne <- data.frame(saver@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
alra.tsne <- data.frame(alra@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
library(ggplot2)
library(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my.cols <- gg_color_hue(length(unique(ident2)))


xshift <- 14
p1 <- ggplot(obs.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), size = 0.5) + guides(colour = FALSE) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 20, face = "plain")) + annotate(x=max(obs.tsne$tSNE_1)-xshift,y=min(obs.tsne$tSNE_2), label=sprintf("OOB Error = %.1f%%", rfr$prediction.error*100),geom="text") + labs(title = "Observed")
p2 <- ggplot(alra.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), size = 0.5)  + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 20, face = "plain")) +scale_colour_manual(name = "Cell type", values = my.cols)  + guides(colour=FALSE)+annotate(x=max(alra.tsne$tSNE_1)-xshift,y=min(alra.tsne$tSNE_2), label=sprintf("OOB Error = %.1f%%", rfr.alra$prediction.error*100),geom="text") + labs(title='ALRA')
p3 <- ggplot(saver.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), size = 0.5)  + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 20, face = "plain"), legend.position="bottom", legend.title=element_blank()) +scale_colour_manual(name = "Cell type", values = my.cols)  + guides(colour=guide_legend(override.aes=list(size=10), nrow=4))+annotate(x=max(saver.tsne$tSNE_1)-xshift,y=min(saver.tsne$tSNE_2), label=sprintf("OOB Error = %.1f%%", rfr.saver$prediction.error*100),geom="text")+ labs(title='SAVER')
g1<-plot_grid(p1,p2,p3+theme(legend.position="none"),nrow=1)
ggsave(g1,filename="figs/Hrvatin-tsne-subset.pdf",width=12,height=4)


gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}


g2 <- gglegend(p2)
ggsave(g2,filename="figs/Hrvatin-tsne_legend.pdf",width=12,height=2)
