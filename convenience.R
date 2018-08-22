
library(grid)
library(ggthemes)
theme_Publication <- function(base_size=14, base_family="Helvetica") {
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(size = rel(0.8)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))

}

pca_tsne <- function(A_norm) {
#center the columns
	A_rm <- colMeans(A_norm)
	A_norm_c <- sweep(A_norm, 2,A_rm) 
	print("Computing SVD")
	library(rsvd)
	fastDecomp <- rsvd(A_norm_c, 50,q=3);
	print("Computing t-SNE")
	PCs<- fastDecomp$u %*% diag(fastDecomp$d);
	setwd('/data/george/Research_Local/FIt-SNE/')
	source('fast_tsne.R')
	tsne.out <- fftRtsne(as.matrix(PCs[,1:50]),data_path = "temp/data.dat", result_path = "temp/result.dat",rand_seed = 4,ann_not_vptree=FALSE)

	list(tsne=tsne.out,PCs=PCs)
}
