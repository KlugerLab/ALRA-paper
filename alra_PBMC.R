##=======================================================##
## Purified PBMC bulk RNA-seq
## ALRA
##=======================================================##



# Load data ---------------------------------------------------------------


library(tidyverse)
library(Matrix)
library(cowplot)
library(gridExtra)
library(rsvd)
library(SAVER)
library(Rmagic)
library(scImpute)
source('ClstAvgExp_facetted.R')
library(reshape2)
library(ggplot2)

source("../ALRA/alra.R")


## Load scRNA-seq data

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
  dir.create('data/purified_pbmcs', showWarnings =F)
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
    genenames <- read.delim(sprintf('%s/filtered_matrices_mex/hg19/genes.tsv', fn2), sep = '\t',header = FALSE)[,2]
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



# Load bulk ---------------------------------------------------------------




## Bulk RNA-seq data

# read in data
fn <-  'data/GSE64655_Normalized_transcript_expression_in_human_immune_cells.txt.gz'
if ( !file.exists(fn)) {
    download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64655&format=file&file=GSE64655%5FNormalized%5Ftranscript%5Fexpression%5Fin%5Fhuman%5Fimmune%5Fcells%2Etxt%2Egz',fn)
}

PBMC_bulk <- read.table( "data/GSE64655_Normalized_transcript_expression_in_human_immune_cells.txt.gz", sep = "\t", header = T, stringsAsFactors = F, skip = 4, quote="\"")

rownames(PBMC_bulk) <- PBMC_bulk$Gene.ID


# scRNA-seq and bulk genes
mean(colnames(A) %in% PBMC_bulk$Gene.Symbol)

# read in scRNA-seq genes and change to ensemble ID
scRNA_genes <- read.table(
  "data/purified_pbmcs/b_cells_unzipped/filtered_matrices_mex/hg19/genes.tsv",
  sep = "\t", header = F, stringsAsFactors = F
)
rownames(scRNA_genes) <- scRNA_genes[,1]
all(scRNA_genes[,2] == colnames(A))

colnames(A) <- scRNA_genes[,1]
mean(colnames(A) %in% PBMC_bulk$Gene.ID)


# Normalize and filter data -----------------------------------------------


num_of_genes_in_cell <- rowSums(A>0)
num_of_cells_in_gene <- colSums(A>0)
keep_cells = which(num_of_genes_in_cell > 400)
keep_genes = which(num_of_cells_in_gene > 100)

A_[1:5,1:5]
A_ <- A[keep_cells,keep_genes]
A__genenames <- scRNA_genes
mean(colnames(A_) %in% PBMC_bulk$Gene.ID)

labels_ <- labels[keep_cells]
A_norm <- normalize_data(A_)

rm(A) # Clear some memory
A_ <- Matrix(A_, sparse = T)
A_norm <- Matrix(A_norm, sparse = T)


# write A_ into file for scImpute and DCA
fn =  "data/purified_PBMC_filtered.csv"
if ( !file.exists(fn)) {
write.table( as.matrix(t(A_)),fn, sep = ",", quote = F, row.names = T,
            col.names = T
)
}else {
    print("Using existing purified PBMC filtered CSV instead of writing new one.")
}



# merge T cell labels
scRNA_labels <- labels_
scRNA_labels[labels_ %in% c("cd4_helper","regulatory_t","naive_t","memory_t","cytotoxic_t","naive_cytotoxic")] <- "t_cells"
unique(scRNA_labels)


# match bulk label & scRNA label, in SAME order
bulklabel2do <- list(c("HD30_B_d0","HD31_B_d0"),
                     c("HD30_Mono_d0","HD31_Mono_d0"),
                     c("HD30_T_d0","HD31_T_d0"),
                     c("HD30_NK_d0","HD31_NK_d0"))
sclabel2do <- c("b_cells","cd14_monocytes","t_cells","cd56_nk")
n_label2do <- length(sclabel2do)


# get zero genes for each label
getZeroGenes <- function(bulk.data, sc.data){
  common.genes <- intersect(rownames(bulk.data), colnames(sc.data))
  bulk.data <- bulk.data[common.genes,]
  zero.genes <- rownames(bulk.data)[rowSums(bulk.data) == 0]
  cat("# zero genes: ", length(zero.genes), "\n", sep = "")
  cat("# zero entries: ", sum(sc.data[,zero.genes] == 0), "\n", sep = "")
  return(zero.genes)
}

PBMC_bulk_zero_genes <- list()
for(i in 1:n_label2do){
  PBMC_bulk_zero_genes[[i]] <- getZeroGenes(
    bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
    sc.data = A_norm[scRNA_labels == sclabel2do[i],]
  )
}
names(PBMC_bulk_zero_genes) <- sclabel2do

length(unlist(PBMC_bulk_zero_genes))
length(unique(unlist(PBMC_bulk_zero_genes)))

for(i in 1:(n_label2do-1)){
  for(j in (i+1):n_label2do){
    cat(
      sclabel2do[i], " ", sclabel2do[j], ": ", 
      length(intersect(PBMC_bulk_zero_genes[[i]], PBMC_bulk_zero_genes[[j]])),
      "zero genes in common\n", sep = ""
    )
  }
}


# check zero ratio originally
for(i in 1:n_label2do){
  cat(
    mean(A_[scRNA_labels == sclabel2do[i],PBMC_bulk_zero_genes[[i]]] == 0), 
    "\n", sep = ""
  )
}



# check zero gene names
lapply(PBMC_bulk_zero_genes[2], FUN = function(x,info){
  sapply(x, FUN = function(xx,info){
    info[info[,1] == xx,2]
  },info = info)
}, info = scRNA_genes)


# write zero genes in table
write.table(
  PBMC_bulk[unique(unlist(PBMC_bulk_zero_genes)), c("Gene.ID", "Gene.Symbol", unlist(bulklabel2do))],
  file = "data/PBMC_bulk_day0_zero_genes2.tsv", sep = "\t", row.names = F, col.names = T, quote = F
)



# Completion by various methods-----------------------------------------------

# ALRA
set.seed(3)
k_choice <- choose_k(A_norm, q=10)
system.time(alra.result <- alra(as.matrix(A_norm),k_choice$k))
A_norm_alra <- alra.result[[3]]

diffs <- k_choice$d[1:(length(k_choice$d)-1)] - k_choice$d[2:length(k_choice$d)] 
df <- data.frame(x=1:99,y=diffs)[10:99,]
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=2)  + geom_line()+ geom_vline(xintercept=k_choice$k)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10))
g

ggsave('figs/spectrum_purified_pbmcs.pdf', width=4,height=2, g)
print(sprintf('Chose k=%d',k_choice$k))

# SAVER on subsample data with all genes
set.seed(3)
fn <- 'data/saver_result_10000cells_seed3.RData'
if ( !file.exists(fn)) {
        print(sprintf("Generating %s\n", fn));
        subsample_saver <- sample(nrow(A_), size = 10000, replace = F)
        A_saver <- saver(t(A_[subsample_saver,]))
        save(list = c("subsample_saver","A_saver"), file =fn)
}else{
    print(sprintf("Loading %s\n", fn));
    load(fn)
}
A_saver_sampled <- t(sample.saver(A_saver, seed=0))


# SAVER on all cells with zero genes only for calculating zero preservation
# ratio
fn =  "data/saver_result_433genes.RData"
if ( !file.exists(fn)) {
    print(sprintf("Generating %s", fn))
    set.seed(3)
    A_saver_zero_genes <- saver(
      t(A_),
      pred.genes = which(colnames(A_) %in% 
                         unique(unlist(PBMC_bulk_zero_genes))), 
      pred.genes.only = T)
    save(list = c("A_saver_zero_genes"), file =fn)
}else{
    print(sprintf("Loading %s", fn));
    load(fn)
}
A_saver_sampled_zero_genes <- t(sample.saver(A_saver_zero_genes, seed=3))



fn ='data/purified_pbmc_saver.RDS';
if ( !file.exists(fn)) {
    print(sprintf("Generating %s\n", fn));
    system.time(A_saver_markers <- saver(t(A_),pred.genes = which( scRNA_genes[keep_genes,2] %in% c("CD19","PAX5", "NCAM1","CD4", "CD8A")),pred.genes.only = TRUE))
    saveRDS(A_saver_markers,fn)
}else{
    print(sprintf("Loading %s\n", fn));
    A_saver_markers <- readRDS(fn)
}
A_saver_sampled_markers <- t(sample.saver(A_saver_markers,seed=3))
colnames(A_saver_sampled_markers) <- scRNA_genes[colnames(A_saver_sampled_markers),2]
A_saver_sampled_markers <- log(A_saver_sampled_markers + 1)

fn ='data/purified_pbmc_magic.RDS';
if ( !file.exists(fn)) {
        print(sprintf("Generating %s\n", fn));
	system.time(	magic_result <- magic(A_norm))
	saveRDS(magic_result,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        magic_result <- readRDS(fn)
}
magic_norm <- magic_result$result
colnames(magic_norm) <- colnames(A_norm)


# scImpute
set.seed(3)
fno =  "data/pbmc_scimpute_out/scimpute_count.csv"
if ( !file.exists(fno)) {
    print(sprintf("Generating %s\n", fno));
        fni <- "data/purified_PBMC_filtered_10000.csv"
        write.table( as.matrix(t(A_[subsample_saver,])), 
                    fni, sep = ",", quote = F, row.names = T, col.names = T)
        scimpute_result <- scimpute(
          count_path =fni, 
          out_dir = "data/pbmc_scimpute_out/", 
          Kcluster = 10
          )
}
scimpute_norm <- as.matrix(read.table( fno, sep = ",", header = T, row.names = 1, stringsAsFactors = F))

mean(scimpute_norm == 0)

scimpute_norm <- Matrix(scimpute_norm, sparse = T)
scimpute_norm <- t(normalize_data(t(scimpute_norm)))

aggregate(data.frame(scimpute_norm['ENSG00000153563',]), by=list(scRNA_labels[subsample_saver]), FUN=function(x) mean (x>0))
aggregate(data.frame(scimpute_norm['ENSG00000010610',]), by=list(labels_[subsample_saver]), FUN=function(x) mean (x>0))


# DCA
fno =  "data/pbmc_dca_out/mean.tsv"
if ( !file.exists(fno)) {
    print(sprintf("Generating %s", fno))
    system('dca data/purified_PBMC_filtered.csv data/pbmc_dca_out')
}
dca_norm <- as.matrix(read.table(fno, sep = "\t", header = T, row.names = 1, stringsAsFactors = F))
dca_norm <- t(normalize_data(t(dca_norm)))

# Generate Dot Plot ---------------------------------------------------------------

dca_norm_t <- t(dca_norm)

# create dataset
cois <- rev(c("b_cells","cd56_nk","cd4_helper","cytotoxic_t"))
subset_types <- labels_ %in% cois
gois <- c("PAX5", "NCAM1", "CD4", "CD8A")
idois <- sapply(gois, function(x,info){
  info[info[,2] == x,1]
}, info = scRNA_genes)

# Because scImpute uses less genes than the rest, we need to add it
# separately. So, let's first make the dotplot for ALRA, MAGIC, and SAVER.
# Then we'll add on scImpute and DCA
toplot <- data.frame(
  orig=as.matrix(A_norm)[subset_types,idois],
  alra=as.matrix(A_norm_alra)[subset_types,idois],
  magic=magic_norm[subset_types,idois],
  saver=A_saver_sampled_markers[subset_types,gois]
)


# change colnames
for(i in 1:length(idois)){
  colnames(toplot) <- gsub(idois[i], gois[i], colnames(toplot))
}
# change cell type names
some_levels_renamed <- factor(labels_[subset_types], levels = cois)
levels(some_levels_renamed) <- rev(c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells"))
g.data1 <- ClstAvgExp_facetted(
  t(toplot), some_levels_renamed, levels(some_levels_renamed), gois, 
  c("orig","alra","magic","saver"), geneGroup = rep(gois, each = 4), return.data = T
)
str(labels_[subsample_saver][labels_[subsample_saver]%in% cois])
# add scImpute result in dot plot
some_levels_renamed_subsample <- factor(
  labels_[subsample_saver][labels_[subsample_saver] %in% cois], levels = cois
)
levels(some_levels_renamed_subsample) <- rev(c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells"))
toplot <- data.frame(
  scimpute = as.matrix(t(scimpute_norm)[labels_[subsample_saver] %in% cois,idois])
)
for(i in 1:length(idois)){
  colnames(toplot) <- gsub(idois[i], gois[i], colnames(toplot))
}
g.data2 <- ClstAvgExp_facetted(
  t(toplot), some_levels_renamed_subsample, levels(some_levels_renamed_subsample), gois, c("scimpute"), geneGroup = rep(gois, each = 1), return.data = T
)
g.data2

# add dca result in dot plot
toplot <- data.frame(dca=dca_norm_t[subset_types,idois])
for(i in 1:length(idois)){
  colnames(toplot) <- gsub(idois[i], gois[i], colnames(toplot))
}
g.data3 <- ClstAvgExp_facetted(
  t(toplot), some_levels_renamed, levels(some_levels_renamed), gois, 
  c("dca"), geneGroup = rep(gois, each = 1),  return.data = T
)


# Combine and plot
g.data <- rbind(g.data1, g.data2, g.data3)
g <- ggplot(data = g.data, aes(x = method, y = (Cluster), size = NumCells, col = AverageExpression)) + 
  geom_point() + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(.~genegroup,scales="free_x") +
  ylab("Cluster") + scale_colour_gradient(low = "blue", high = "red", guide_colourbar(title = "Average\nExpression\n")) + 
  scale_size(name = "% Cells > 0", range = c(0, 13)) + 
  guides(size = guide_legend(order=1), color = guide_colourbar(order = 2)) + 
  theme(axis.text = element_text(size = 19), axis.title = element_text(size = 17), 
        axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size=19),
  legend.text = element_text(size=14),
  legend.title = element_text(size=19),
  )
g
ggsave('figs/PBMC_dotplot.pdf', plot = g, width = 15, height = 6)



# Generate Density plots ---------------------------------------------------------------



goi <- gois[1]
idoi <- idois[goi]
some_levels_renamed <- factor(some_levels_renamed, levels = c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells"))
glist <- lapply( idois, function (goi)  
{
  x <- alra.result[[1]][subset_types,goi]
  thresh <- quantile(x,0.001)
  ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
    geom_vline(xintercept = c(thresh, -thresh), linetype="dotted", size=1, color =c("blue","red")) +
    scale_x_continuous(breaks=c(0))+
    theme(axis.line.y=element_blank())
})
(g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1)))
ggsave('figs/density_plots1.pdf', width=8,height=1.5, g)

legend <- get_legend(glist[[1]])
ggsave('figs/pbmc_density_legend.pdf', width=5,height=0.8, legend)

glist <- lapply( idois, function (goi)  
{
  x <- magic_result$result[subset_types,goi]
  ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
    #scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
    theme(axis.line.y=element_blank())
})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))

ggsave('figs/density_plots_magic.pdf', width=8,height=1.5, g)

glist <- lapply( gois, function (goi)  
{
  #nonzero <- A_saver_sampled[subset_types,goi] >0
  #x <- A_saver_sampled[subset_types,goi]
    x<-A_saver_sampled_markers[subset_types,goi]
  ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
    geom_vline(xintercept = c(min(x), -min(x)), linetype="dotted", size=1) +
    #scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
    theme(axis.line.y=element_blank())
})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))

ggsave('figs/density_plots_saver.pdf', width=8,height=1.5, g)


# alra
glist <- lapply(idois, function (goi)  
{
  x <- alra.result[[1]][subset_types,goi]
  ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(
      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), 
      legend.title=element_blank(),legend.position="bottom"
    ) +
    geom_vline(xintercept = c(min(x), -min(x)), linetype="dotted", size=1) +
    scale_x_continuous(breaks=c(0))+
    theme(axis.line.y=element_blank())
})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))

ggsave("figs/PBMC_density_alra.pdf", plot = g, width=8, height=1.5)


# scImpute
some_levels_renamed_subsample <- factor(
  some_levels_renamed_subsample, c("B Cells", "NK Cells","Helper T Cells","Cytotoxic T Cells")
)
length(some_levels_renamed_subsample)

#head(some_levels_renamed_subsample)
#head(labels_[subsample_saver][labels_[subsample_saver] %in% cois])
glist <- lapply(idois, function (goi)  
{
  x <- t(scimpute_norm)[labels_[subsample_saver] %in% cois,goi]
  ggplot() + aes(x, fill=some_levels_renamed_subsample) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(
      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), 
      legend.title=element_blank(),legend.position="bottom"
    ) + 
    #scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
    theme(axis.line.y=element_blank())
})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))
ggsave('figs/PBMC_density_scImpute.pdf', plot = g, width=8,height=1.5)


# DCA
glist <- lapply(idois, function (goi)  
{
  x <- dca_norm_t[subset_types,goi]
  ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
    theme(
      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), 
      legend.title=element_blank(),legend.position="bottom"
    ) + 
  #  scale_x_continuous(breaks=c(0), limits = c(min(x,0),quantile(x,0.99))) + 
    theme(axis.line.y=element_blank())
})
g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1))

ggsave("figs/PBMC_density_DCA.pdf", plot = g, width=8, height=1.5)





# Define biological zeros based on bulk RNA-seq -----------------------------------------------


## Compare with bulk

# function to evaluate biological zero
evalBioZero <- function(bulk.data, sc.data, impute.data, return.genes = F){
  common.genes <- intersect(rownames(bulk.data), colnames(sc.data))
  bulk.data <- bulk.data[common.genes,]
  zero.genes <- rownames(bulk.data)[rowSums(bulk.data) == 0]
  cat("# zero genes: ", length(zero.genes), "\n", sep = "") 
  num.zeros <- length(zero.genes) * nrow(sc.data)
  sc.zero.entry <- sum(sc.data[,zero.genes] == 0) / num.zeros 
  sc.zero.gene <- sum(colSums(sc.data[,zero.genes]) == 0) / length(zero.genes)
  # cat("Nonzero ratio of zero genes in scRNA: ", sc.nonzero, "\n", sep = "") 
  res.entry <- sum(sc.data[,zero.genes] == 0 & impute.data[,zero.genes] <= 0) / sum(sc.data[,zero.genes] == 0)
  res.gene <- sum(colSums(impute.data[,zero.genes] > 0) == 0 & colSums(sc.data[,zero.genes]) == 0) / 
    sum(colSums(sc.data[,zero.genes]) == 0) 
  if(return.genes){
    return(list(
      zero.genes = zero.genes,
      res = c(sc.zero.entry, sc.zero.gene, res.entry, res.gene)
    ))
  } else {
    return(c(sc.zero.entry, sc.zero.gene, res.entry, res.gene))
  }
}


# evaluate each cell type
# the 3rd column represents zero preserve ratio
zero_preserve <- matrix(0, nrow = n_label2do, ncol = 4)
rownames(zero_preserve) <- sclabel2do
for(i in 1:n_label2do){
  zero_preserve[i,] <- evalBioZero(
    bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
    sc.data = A_norm[scRNA_labels == sclabel2do[i],],
    impute.data = A_norm_alra[scRNA_labels == sclabel2do[i],]
  )
}

zero_complete_alra <- rep(0, n_label2do)
names(zero_complete_alra) <- sclabel2do
for(i in 1:n_label2do){
  zero_complete_alra[i] <- mean(A_norm[scRNA_labels == sclabel2do[i],] == 0 & A_norm_alra[scRNA_labels == sclabel2do[i],] > 0)
}

# MAGIC
zero_preserve_magic <- matrix(0, nrow = n_label2do, ncol = 4)
rownames(zero_preserve_magic) <- sclabel2do
for(i in 1:n_label2do){
  zero_preserve_magic[i,] <- evalBioZero(
    bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
    sc.data = A_norm[scRNA_labels == sclabel2do[i],],
    impute.data = magic_norm[scRNA_labels == sclabel2do[i],]
  )
}
zero_complete_magic <- zero_complete_alra
for(i in 1:n_label2do){
  zero_complete_magic[i] <- mean(A_norm[scRNA_labels == sclabel2do[i],] == 0 & magic_norm[scRNA_labels == sclabel2do[i],] > 0)
}

# SAVER
zero_preserve_saver <- matrix(0, nrow = n_label2do, ncol = 4)
rownames(zero_preserve_saver) <- sclabel2do
for(i in 1:n_label2do){
  zero_preserve_saver[i,] <- evalBioZero(
    bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
    sc.data = A_norm[scRNA_labels == sclabel2do[i],],
    impute.data = A_saver_sampled_zero_genes[scRNA_labels == sclabel2do[i],]
  )
}
zero_complete_saver <- zero_complete_alra
for(i in 1:n_label2do){
  zero_complete_saver[i] <- mean((A_norm[subsample_saver,][scRNA_labels[subsample_saver] == sclabel2do[i],colnames(A_saver_sampled)] == 0) & (A_saver_sampled[scRNA_labels[subsample_saver] == sclabel2do[i],] > 0))
}

# scImpute
zero_preserve_scimpute_sub <- rep(0,n_label2do)
names(zero_preserve_scimpute_sub) <- sclabel2do
for(i in 1:n_label2do){
  zero_preserve_scimpute_sub[i] <- 
    sum(scimpute_norm[PBMC_bulk_zero_genes[[i]],scRNA_labels[subsample_saver] == sclabel2do[i]] == 0) / 
    sum(A_[subsample_saver,][scRNA_labels[subsample_saver] == sclabel2do[i],PBMC_bulk_zero_genes[[i]]] == 0)
}
zero_complete_scimpute <- zero_complete_alra
for(i in 1:n_label2do){
  #zero_complete_scimpute[i] <- zero_preserve_scimpute_sub[i] - 
  #  mean(scimpute_norm[,scRNA_labels[subsample_saver] == sclabel2do[i]] == 0)
  zero_complete_scimpute[i] <- mean(
    (A_norm[subsample_saver,] [scRNA_labels[subsample_saver] == sclabel2do[i],]) & (t(scimpute_norm)[scRNA_labels[subsample_saver] == sclabel2do[i],] > 0)
  ) 
}

for(i in 1:4){
    cat(sprintf('%.2f, %.2f\n', zero_preserve[i,3], zero_complete_alra[i]))
}
for(i in 1:4){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_magic[i,3], zero_complete_magic[i]))
}
for(i in 1:4){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_saver[i,3], zero_complete_saver[i]))
}
for(i in 1:4){
    cat(sprintf('%.2f, %.2f\n', zero_preserve_scimpute_sub[i], zero_complete_scimpute[i]))
}




################################
# Gating 
################################
A_norm_alra_2 <- A_norm_alra
colnames(A_norm_alra_2)<- scRNA_genes[keep_genes,2]
celltypes  <- list(
                  # myeloidDC = list(
                  #              pos_marks = c('ITGAX', 'HLA-DRA'),
                  #              neg_marks = c('CD3G', 'CD19', 'MS4A1', 'NCAM1' , 'FCGR3A' )
                  #              ),
                   plasmacytoidDC = list(
                                pos_marks = c('IL3RA', 'HLA-DRA'),
                                neg_marks = c('CD3G', 'CD19', 'MS4A1', 'NCAM1' , 'ITGAX', 'CD14')
                                ),
                   Treg = list(
                                pos_marks = c( 'CD3G', 'CD4', 'IL2RA', 'CCR4', 'FOXP3'),
                                neg_marks = c()
                                ),
                   memoryCD4 = list(
                                pos_marks = c( 'CD3G', 'CD4', 'CCR7', 'SELL'),
                                neg_marks = c( 'CD8A', 'CD8A')
                                ),
                   memoryCD8 = list(
                                pos_marks = c( 'CD3G', 'CD8A', 'CCR7', 'SELL'),
                                neg_marks = c( 'CD4', 'CD4')
                                ),
                   memoryB = list(
                                pos_marks = c( 'CD27', 'CD19'),
                                neg_marks = c( 'CD14', 'CD3G')
                                ),
                   naiveB = list(
                                pos_marks = c('CD19', 'CD19'),
                                neg_marks = c( 'CD14', 'CD27', 'CD3G')
                                ),
                   nkCells = list(
                                pos_marks = c('NCAM1', 'NCAM1'),
                                neg_marks = c('CD3G', 'CD19')
                                ),
                   classical_monocyte = list(
                                pos_marks = c('CD14', 'CD14'),
                                neg_marks = c('CD3G', 'CD19', 'MS4A1', 'NCAM1', 'FCGR3A')
                                ),
                  # nonclassical_monocyte = list(
                  #              pos_marks = c( 'FCGR3A',  'FCGR3A'),
                  #              neg_marks = c('CD3G', 'CD19', 'MS4A1', 'NCAM1')
                  #              ),
                   intermediate_monocyte = list(
                                pos_marks = c('CD14', 'FCGR3A'),
                                neg_marks = c('CD3G', 'CD19', 'MS4A1', 'NCAM1')
                                )
                   )
cell_labels <- rep('none', nrow(A_norm))
for ( celltype in names(celltypes)) {
    pos_marks <- rowSums(A_norm_alra_2[, celltypes[[celltype]]$pos_marks] > 0  ) == length(celltypes[[celltype]]$pos_marks)
    if ( length(celltypes[[celltype]]$neg_marks) >0 ){
        neg_marks <- rowSums(A_norm_alra_2[, celltypes[[celltype]]$neg_marks] > 0  ) ==0
        identified <- pos_marks & neg_marks
    }else{
        identified <- pos_marks
    }
    cell_labels[identified] <- celltype
    print(celltype)
    print(sum(identified))
}

dim(A_)
A_df  <- as.data.frame(as.matrix(A_))
colnames(A_df)  <- scRNA_genes[keep_genes,2]
A_df$gating_labels <- cell_labels
A_df$true_labels <- labels[keep_cells]
length(A_df$true_labels)
length(A_df$gating_labels)



# Get the most correlated cell type
blood_cell_gene_data <- read_tsv('data/rna_blood_cell.tsv')
unique( blood_cell_gene_data$`Blood cell`)
blood_cell_gene_data <- blood_cell_gene_data %>% rename (gene_name = `Gene name`, blood_cell = `Blood cell`) %>% select(gene_name,blood_cell, TPM) %>% pivot_wider(id_cols = blood_cell, names_from =gene_name, values_from =TPM, values_fn = list(TPM = mean   ))
blood_cell_gene_data <- blood_cell_gene_data %>% mutate(blood.cell.type.chosen = case_when( blood_cell == 'classical monocyte' ~ 'classical_monocyte',
                                                                   blood_cell  == 'myeloid DC' ~'myeloidDC',
                                                                   blood_cell  == 'plasmacytoid DC' ~'plasmacytoidDC',
                                                                   blood_cell == 'T-reg' ~ 'Treg',
                                                                   blood_cell == 'memory CD4 T-cell' ~ 'memoryCD4',
                                                                   blood_cell == 'memory CD8 T-cell' ~ 'memoryCD8',
                                                                   blood_cell == 'memory B-cell' ~ 'memoryB',
                                                                   blood_cell == 'naive B-cell' ~ 'naiveB',
                                                                   blood_cell == 'NK-cell' ~ 'nkCells',
                                                                   blood_cell == 'classical monocyte' ~ 'classical_monocyte',
                                                                   blood_cell == 'intermediate monocyte' ~ 'intermediate_monocyte',
                                                                   blood_cell == 'non-classical monocyte' ~ 'nonclassical_monocyte',
                                                                   T ~ NA_character_))  %>%
                                                filter( !is.na(blood.cell.type.chosen))

# Identify genes in common
common_gn <- intersect(colnames(A_df) , colnames(blood_cell_gene_data))
common_types <- intersect( A_df$cell_labels, (blood_cell_gene_data$blood.cell.type.chosen))
A_df_corrs <- cor( A_df[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t, method = 'spearman')
A_df$correlation_labels <- blood_cell_gene_data$blood.cell.type.chosen[ apply(A_df_corrs, 1, which.max)]


A_original_temp <- as.data.frame(as.matrix(A)) %>% setNames(make.names(names(.), unique = TRUE)) 
A_original_temp %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) sum(r>0))
A_original_temp$true_labels  <- labels
A_original_temp %>% group_by(true_labels) %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) mean(r))



table( A_df$true_labels ,A_df$gating_labels, useNA="ifany")
table( A_df$true_labels ,A_df$correlation_labels, useNA="ifany")





################################
# Test scale factor
################################

normalize_data <- function (A, s.f = 10E3) {
  #  Simple convenience function to library and log normalize a matrix
  
  totalUMIPerCell <- rowSums(A);
  if (any(totalUMIPerCell == 0)) {
    toRemove <- which(totalUMIPerCell == 0)
    A <- A[-toRemove,]
    totalUMIPerCell <- totalUMIPerCell[-toRemove]
    cat(sprintf("Removed %d cells which did not express any genes\n", length(toRemove)))
  }
  
  A_norm <- sweep(A, 1, totalUMIPerCell, '/');
  A_norm <- A_norm * s.f
  A_norm <- log(A_norm +1);
}

for(sf in c(1000,5000,10000,20000)){
  A_norm <- normalize_data(A_, s.f = sf)
  A_norm <- Matrix(A_norm, sparse = T)
  set.seed(3)
  k_choice <- choose_k(A_norm, q=10)
  alra.result <- alra(as.matrix(A_norm),k_choice$k)
  A_norm_alra <- alra.result[[3]]
  
  zero_preserve <- matrix(0, nrow = n_label2do, ncol = 4)
  rownames(zero_preserve) <- sclabel2do
  for(i in 1:n_label2do){
    zero_preserve[i,] <- evalBioZero(
      bulk.data = PBMC_bulk[,bulklabel2do[[i]]],
      sc.data = A_norm[scRNA_labels == sclabel2do[i],],
      impute.data = A_norm_alra[scRNA_labels == sclabel2do[i],]
    )
  }
  zero_complete_alra <- rep(0, n_label2do)
  names(zero_complete_alra) <- sclabel2do
  for(i in 1:n_label2do){
    zero_complete_alra[i] <- mean(A_norm[scRNA_labels == sclabel2do[i],] == 0 & A_norm_alra[scRNA_labels == sclabel2do[i],] > 0)
  }
  for(i in 1:4){
    cat(sprintf('%.2f, %.2f\n', zero_preserve[i,3], zero_complete_alra[i]))
  }
  for(i in 1:4){
    cat(sprintf('%.2f\n', zero_preserve[i,3]))
  }
  for(i in 1:4){
    cat(sprintf('%.2f\n', zero_complete_alra[i]))
  }
  
  glist <- lapply( idois, function (goi)
  {
    x <- alra.result[[1]][subset_types,goi]
    thresh <- quantile(x,0.001)
    ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") + 
      geom_vline(xintercept = c(thresh, -thresh), linetype="dotted", size=1, color =c("blue","red")) +
      theme(axis.line.y=element_blank())
  })
  (g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1)))
  ggsave(paste0('figs/density_lra_',sf,'.pdf'), width=8,height=1.5, g)
  
  glist <- lapply( idois, function (goi)
  {
    x <- alra.result[[3]][subset_types,goi]
    ggplot() + aes(x, fill=some_levels_renamed) + geom_density(aes(y=..scaled..), alpha=1/2) + 
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(),legend.position="bottom") +
      #scale_x_continuous(breaks=c(0), limits = c(min(x),quantile(x,0.99)))+
      theme(axis.line.y=element_blank())
  })
  (g <- do.call("grid.arrange", c( lapply(glist, '+', theme(legend.position="none")), nrow=1)))
  ggsave(paste0('figs/density_alra_',sf,'.pdf'), width=8,height=1.5, g)
}





################################
# Correlation with nonzero
################################

nonzero_idx <- which(A_norm != 0)
cor(A_norm[nonzero_idx], A_norm_alra[nonzero_idx], method = "pearson")
cor(A_norm[nonzero_idx], A_norm_alra[nonzero_idx], method = "spearman")
cor(A_norm[nonzero_idx], alra.result[[2]][nonzero_idx], method = "pearson")
cor(A_norm[nonzero_idx], alra.result[[2]][nonzero_idx], method = "spearman")

magic_norm <- as.matrix(magic_norm)
cor(A_norm[nonzero_idx], magic_norm[nonzero_idx], method = "pearson")
cor(A_norm[nonzero_idx], magic_norm[nonzero_idx], method = "spearman")

dca_norm <- t(dca_norm)
cor(A_[nonzero_idx], dca_norm[nonzero_idx], method = "pearson")
cor(A_[nonzero_idx], dca_norm[nonzero_idx], method = "spearman")
