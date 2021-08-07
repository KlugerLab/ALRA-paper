
library(grid) 

theme_tsne <- theme_cowplot() + theme(axis.title = element_blank(), axis.text = element_blank(),
                    axis.ticks = element_blank(), plot.title =
                        element_text(size = 10, face = "plain"),
                    legend.position="none",
                    plot.margin=unit(c(0.05,0.05,0.05,0.05 ), "cm"),
                    axis.line.y.left=element_line(colour="black",size=0.2),
                    axis.line.x.bottom=element_line(colour="black",size=0.2))

#theme_Publication <- function(base_size=14, base_family="Helvetica") {
#      (theme_foundation(base_size=base_size, base_family=base_family)
#       + theme(plot.title = element_text(
#                                         size = rel(1.2), hjust = 0.5),
#               text = element_text(),
#               panel.background = element_rect(colour = NA),
#               plot.background = element_rect(colour = NA),
#               panel.border = element_rect(colour = NA),
#               axis.title = element_text(size = rel(0.8)),
#               axis.title.y = element_text(angle=90,vjust =2),
#               axis.title.x = element_text(vjust = -0.2),
#               axis.text = element_text(),
#               axis.line = element_line(colour="black"),
#               axis.ticks = element_line(),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               legend.key = element_rect(colour = NA),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.key.size= unit(0.2, "cm"),
#               legend.title = element_text(face="italic"),
#               plot.margin=unit(c(10,5,5,5),"mm"),
#               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#               strip.text = element_text(face="bold")
#          ))
#
#}

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




if_not_exist_do <- function( fn, code, codeelse={ cat(sprintf("%s does exist\n", fn)) } ) {
    if ( !file.exists(fn)) {
        print(sprintf("%s does not exist", fn))
          code;
    }else {
        codeelse;
    }
}


annotation_compass <- function(label,
                               position = c('N','NE','E','SE','S','SW','W','NW'),
                               padding = grid::unit(c(0.0,0.25),"line"),  fontsize=9){
  position <- match.arg(position)
  x <- switch (position,
    N = 0.5,
    NE = 1,
    E = 1,
    SE = 1,
    S = 0.5,
    SW = 0,
    W = 0,
    NW = 0
  )
  y <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0,
               SW = 0,
               W = 0.5,
               NW = 1
  )
  hjust <- switch (position,
               N = 0.5,
               NE = 1,
               E = 1,
               SE = 1,
               S = 0.5,
               SW = 0,
               W = 0,
               NW = 0
  )

  vjust <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0,
               SW = 0,
               W = 0.5,
               NW = 1
  )
  f1 <- switch (position,
                   N = 0,
                   NE = -1,
                   E = -1,
                   SE = -1,
                   S = 0,
                   SW = 1,
                   W = 1,
                   NW = 1
  )
  f2 <- switch (position,
                   N = -1,
                   NE = -1,
                   E = 0,
                   SE = 1,
                   S = 1,
                   SW = 1,
                   W = 0,
                   NW = -1
  )
  annotation_custom(grid::textGrob(label,
                                   x=grid::unit(x,"npc") + f1*padding[1] ,
                                   y=grid::unit(y,"npc") + f2*padding[2],
                                   hjust=hjust,vjust=vjust, gp=gpar(fontsize=fontsize)))
}
