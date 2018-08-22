library(Matrix)
library(cowplot)
library(gridExtra)
library(rsvd)
library(SAVER)
source('gene_cluster_heatmap.R')
source('../ALRA/alra.R')
library(reshape2)
library(ggplot2)

if ( !file.exists('data/purified_pbmc.RData')) {
	urls <- list(b_cells='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz', 
		     cd14_monocytes='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_filtered_gene_bc_matrices.tar.gz',
		     cd34='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd34/cd34_filtered_gene_bc_matrices.tar.gz',
		     cd4_helper='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/cd4_t_helper_filtered_gene_bc_matrices.tar.gz',
		     regulatory_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/regulatory_t/regulatory_t_filtered_gene_bc_matrices.tar.gz',
		     naive_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_t/naive_t_filtered_gene_bc_matrices.tar.gz',
		     memory_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/memory_t/memory_t_filtered_gene_bc_matrices.tar.gz',
		     cd56_nk='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/cd56_nk_filtered_gene_bc_matrices.tar.gz',
		     cytotoxic_t='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/cytotoxic_t/cytotoxic_t_filtered_gene_bc_matrices.tar.gz',
		     naive_cytotoxic='http://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_cytotoxic/naive_cytotoxic_filtered_gene_bc_matrices.tar.gz'	)


	A <- c()
	labels = c()
	for (i in seq_along(urls)) {
		print(i)
		label <-names(urls)[i]
		fn <- sprintf('data/purified_pbmcs/%s.tar.gz', label)
		download.file(urls[[i]],fn)
		# untar
		fn2 <- sprintf('data/purified_pbmcs/%s_unzipped' ,label)
		untar(fn, exdir=fn2)
		mtx <- as.matrix((readMM(
					 sprintf('%s/filtered_matrices_mex/hg19/matrix.mtx', fn2))))
		genenames <- read.delim(sprintf('%s/filtered_matrices_mex/hg19/genes.tsv', fn2),
					sep = '\t',header = FALSE)[,2]
		rownames(mtx) <- genenames
		if (i>1 && !all(rownames(mtx) == colnames(A))) {
			error('Trying to concatenate a matrix with different genenames')
		}
		A <- rbind(A,t(mtx))
		labels <- c(labels, rep(label,ncol(mtx)))
	}
	save(A,labels,file='data/purified_pbmc.RData')
}else{
	print("Loading");
	load('data/purified_pbmc.RData')
}


num_of_genes_in_cell <- rowSums(A>0)
num_of_cells_in_gene <- colSums(A>0)
keep_cells = which(num_of_genes_in_cell > 400) 
keep_genes = which(num_of_cells_in_gene > 100)

A_ <- A[keep_cells,keep_genes]
labels_ <- as.factor(labels[keep_cells])

A_norm <- normalize_data(A_);

# Completion ------------------------------------------------------------
set.seed(3)
print("Beginning to choose k")
k_choice <- choose_k(A_norm)

df <- data.frame(x=2:100,y=diff(k_choice$d))[10:99,]
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=2)  + geom_line()+ geom_vline(xintercept=k_choice$k+1)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10))
g
ggsave('figs/spectrum_purified_pbmcs.pdf', width=4,height=2, g)
print(sprintf('Chose k=%d',k_choice$k))

system.time(alra.result <- alra(A_norm,k_choice$k))
A_norm_alra <- alra.result[[3]]
#rm(alra.result); gc();

if ( !file.exists('data/purified_pbmc_saver.RDS')) {
	system.time(A_saver <- saver(t(A_),pred.genes = which( colnames(A_) %in% c("CD19","CR2", "NCAM1","CD4", "CD34","CD8A","IL2RA","CD14")),pred.genes.only = TRUE))
	saveRDS(A_saver,'data/purified_pbmc_saver.RDS')
}else{
	A_saver <- readRDS('data/purified_pbmc_saver.RDS')
}

A_saver_sampled <- t(sample.saver(A_saver,seed=3))


if ( !file.exists('data/purified_pbmc_magic.RDS')) {
	library(Rmagic)
	system.time(	magic_result <- magic(A_norm))
	saveRDS(magic_result,'data/purified_pbmc_magic.RDS')
}else{
	magic_result <- readRDS('data/purified_pbmc_magic.RDS')
}


cois <- rev(c("b_cells","cd56_nk","cd4_helper","cytotoxic_t"))
subset_types <- labels_ %in% cois 
gois <- c(    "CR2", "NCAM1",  "CD4","CD8A")
toplot <- data.frame(orig=A_norm[subset_types,gois] , 
		     alra=A_norm_alra[subset_types,gois],
			saver=(A_saver_sampled[subset_types,gois]),
		     magic=magic_result$result[subset_types,gois]
                    )
some_levels_renamed <- factor(droplevels(labels_[subset_types]), levels=cois)
levels(some_levels_renamed) <- rev(c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells"))
ClstAvgExp_facetted(t(toplot),some_levels_renamed,levels(some_levels_renamed),gois, c("orig","alra","saver","magic"), geneGroup = rep(gois, each = 4), scale = 13)
ggsave('figs/purified_pbmc_dotplot_completion.pdf',width=10,height=4.25)


some_levels_renamed <- factor(some_levels_renamed,levels=c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells"))
glist <- lapply( gois, function (goi)  
		{
x <- alra.result[[1]][subset_types,goi]
ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
	geom_vline(xintercept = c(min(x), -min(x)), linetype="dotted", size=1) +
	scale_x_continuous(breaks=c(0))+
	theme(axis.line.y=element_blank())
		})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))
ggsave('figs/density_plots1.pdf', width=8,height=1.5, g)

legend <- get_legend(glist[[1]])
ggsave('figs/density_legend.pdf', width=5,height=0.8, legend)

glist <- lapply( gois, function (goi)  
		{
x <- magic_result$result[subset_types,goi]
ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
	scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
	theme(axis.line.y=element_blank())
		})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))
ggsave('figs/density_plots_magic.pdf', width=8,height=1.5, g)

glist <- lapply( gois, function (goi)  
		{
#nonzero <- A_saver_sampled[subset_types,goi] >0
x <- A_saver_sampled[subset_types,goi]
ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
	geom_vline(xintercept = c(min(x), -min(x)), linetype="dotted", size=1) +
	scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
	theme(axis.line.y=element_blank())
		})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))
ggsave('figs/density_plots_saver.pdf', width=8,height=1.5, g)


#goi <- 'CR2'
#aggregate(A_norm[,goi], by=list(labels_),FUN=function(x) round(c(mean=mean(x), percent=sum(x>0)/length(x)),3))
#aggregate(A_norm_alra[,goi], by=list(labels_[]),FUN=function(x) round(c(mean=mean(x), percent=sum(x>0)/length(x)),3))
#aggregate(A_saver_sampled[,goi], by=list(labels_[]),FUN=function(x) round(c(mean=mean(x), percent=sum(x>0)/length(x)),3))
#aggregate(A_norm_alra[subset_types,goi], by=list(labels_[subset_types]),FUN=function(x) round(c(mean=mean(x), percent=sum(x>0)/length(x)),3))
#aggregate(A_saver_sampled[subset_types,goi], by=list(labels_[subset_types]),FUN=function(x) round(c(mean=mean(x), percent=sum(x>0)/length(x)),3))
#
