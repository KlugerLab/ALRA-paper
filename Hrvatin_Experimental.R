## Hrvatin t-sne analysis with DCA results
# SAVER-data folder can be found here: https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0 
library(ranger)
library(ggrepel)
library(gridExtra)
library(ggplot2) 
library(tidyverse)
library(Seurat) #Seurat Version 2
#devtools::~/Research_Local/seurat/seurat2/seurat
library(ggplot2)
library(cowplot)
library(Rmagic)
source('convenience.R')
source("../ALRA/alra.R")


## Load & filter data

if (!file.exists('data/hrvatin_full.rds')) {
    download.file('https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AACz-gmkCDYojbQkrXmnA1f9a/GSE102827_merged_all_raw.csv.gz?dl=1', 'data/GSE102827_merged_all_raw.csv.gz')
  x <- read.csv("data/GSE102827_merged_all_raw.csv.gz", header = TRUE, row.names = 1, check.names = FALSE)
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
# write filtered data for DCA
    fn <- 'data/GSE102827_filtered_raw.csv'
    write.csv(x, file =fn, row.names = T, col.names = T, quote = F)
}else{
  print("Reading")
  x <- readRDS('data/hrvatin_full.rds')
}




# Get the cell type namevs
ct.dat <- read.csv("data/hrvatin_celltypes.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
ct.dat2$celltype[!is.na(ct.dat2$subtype)] <- 
  ct.dat2$subtype[!is.na(ct.dat2$subtype)]
table(ct.dat2$celltype)
ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]
ident <- ct.dat3[colnames(x), 4]


x.sub <- x[, !is.na(ident)]
ident2 <- ct.dat3[colnames(x.sub), 4]
ident2 <- factor(ident2, levels = unique(ident2)[order(tolower(unique(ident2)))])

x_t <- t(x)
totalUMIPerCell <- rowSums(x_t);
x_norm_sum <- sweep(x_t, 1, totalUMIPerCell, '/');
x_norm_sum <- x_norm_sum * 10E3
x_norm_sum_log <- log(x_norm_sum +1);




seuratify  <- function(xx) {
    xx <- CreateSeuratObject(xx, project = "alra")
    xx <- FindVariableFeatures(xx, selection.method = "mvp")
    xx <- ScaleData(xx)
    return(xx)
}
rangerify  <-  function(xx,xx.ident) {
    df <- data.frame(t(xx@assays$RNA[VariableFeatures(xx),]), Y=as.factor(xx.ident))
    rfr.alra <- ranger(Y ~ ., df)
    rfr.alra$prediction.error
}



################################
# Typical ALRA 
################################
set.seed(3)
system.time(k_choice_sum_log <- choose_k(x_norm_sum_log))
system.time(x_norm_sum_log_alra_result <- alra(x_norm_sum_log,k_choice_sum_log$k)) 
x.lra.sum.log <- as.matrix(t(x_norm_sum_log_alra_result[[1]]))
x.alra.sum.log <- as.matrix(t(x_norm_sum_log_alra_result$A_norm_rank_k_cor_sc))
colnames(x.alra.sum.log) <- colnames(x)
x.alra.sum.log.sub <- x.alra.sum.log[,!is.na(ident)]
set.seed(3)
alra.sum.log.seurat  <-  seuratify( x.alra.sum.log.sub)
alra.sum.log.ranger <-   rangerify(alra.sum.log.seurat, ident2)



system.time(result.completed <- alra(x_norm_sum_log,200)) 
x.alra.sum.log.k200 <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
colnames(x.alra.sum.log.k200) <- colnames(x)
x.alra.sum.log.k200.sub <- x.alra.sum.log.k200[,!is.na(ident)]
set.seed(3)
alra.sum.log.k200.seurat  <-  seuratify( x.alra.sum.log.k200.sub)
alra.sum.log.k200.ranger <-   rangerify(alra.sum.log.k200.seurat, ident2)




################################
# ALRA on Raw 
################################



# Compute ALRA
set.seed(3)
system.time(k_choice_sum <- choose_k(x_norm_sum))
system.time(result.completed <- alra(x_norm_sum,k_choice_sum$k)) 
x.alra.sum <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
colnames(x.alra.sum) <- colnames(x)
x.alra.sum.sub <- log(1+x.alra.sum[,!is.na(ident)])
set.seed(3)
alra.sum.seurat  <-  seuratify(x.alra.sum.sub)
alra.sum.ranger  <-  rangerify(alra.sum.seurat, ident2)





################################
# ALRA on Raw but using k=100 
################################
system.time(result.completed <- alra(x_norm_sum,100)) 
x.alra.k100.sum <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
colnames(x.alra.k100.sum) <- colnames(x)
x.alra.k100.sum.sub <- log(1+x.alra.k100.sum[,!is.na(ident)])
alra.k100.sum.seurat  <-  seuratify(x.alra.k100.sum.sub)
set.seed(4)
alra.sum.k100.ranger  <-  rangerify(alra.k100.sum.seurat, ident2)
#0.044 with seed 3, and with 4



################################
# ALRA on Raw but using k=200 
################################
system.time(result.completed <- alra(x_norm_sum,200)) 
x.alra.k200.sum <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
colnames(x.alra.k200.sum) <- colnames(x)
x.alra.k200.sum.sub <- log(1+x.alra.k200.sum[,!is.na(ident)])
alra.k200.sum.seurat  <-  seuratify(x.alra.k200.sum.sub)
set.seed(3)
alra.sum.k200.ranger  <-  rangerify(alra.k200.sum.seurat, ident2)


################################
# ALRA on Raw but using k=500 
################################
system.time(result.completed <- alra(x_norm_sum,500)) 
x.alra.k500.sum <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
colnames(x.alra.k500.sum) <- colnames(x)
x.alra.k500.sum.sub <- log(1+x.alra.k500.sum[,!is.na(ident)])
alra.k500.sum.seurat  <-  seuratify(x.alra.k500.sum.sub)
set.seed(3)
alra.sum.k500.ranger  <-  rangerify(alra.k500.sum.seurat, ident2)

alra.sum.log.ranger
alra.sum.log.k200.ranger
alra.sum.ranger
alra.sum.k100.ranger
alra.sum.k200.ranger
alra.sum.k500.ranger


#alra.sum.log.seurat <- CreateSeuratObject(x.alra.sum.log.sub, project = "alra")
#alra.sum.log.seurat <- FindVariableFeatures(alra.sum.log.seurat, selection.method = "mvp")
#alra.sum.log.seurat <- ScaleData(alra.sum.log.seurat)
#df <- data.frame(t(x.alra.sum.log.sub[VariableFeatures(alra.sum.log.seurat),]), Y=as.factor(ident2))
#rfr.alra <- ranger(Y ~ ., df)
#rfr.alra$prediction.error

################################
# Explore the eigenvectors 
################################
x.sum.log.svd  <- rsvd( x_norm_sum_log, 300, q = 10)
x.sum.svd  <- rsvd( x_norm_sum, 300, q = 10)

x.sum.log.svd$d[1:20]/x.sum.log.svd$d[1]
x.sum.svd$d[1:20]/x.sum.svd$d[1]
x.sum.svd.v  <- apply(abs(x.sum.svd$v[,1:10]),2,sort, decreasing=T)[1:25,] %>% as.data.frame
x.sum.svd.v.genenames  <- apply(abs(x.sum.svd$v[,1:10]),2, function( x) colnames(x_norm_sum)[order(x, decreasing=T)[1:25]])
x.sum.log.svd.v  <- apply(abs(x.sum.log.svd$v[,1:10]),2,sort, decreasing=T)[1:25,]

library(dplyr)
library(tidyr)
x.sum.log.svd.v.longer$name  <- factor(x.sum.log.svd.v.longer$name, levels = sprintf('V%d', 1:10)) 


order, decreasing=T)[1:25,] %>% as.data.frame
colnames(x_norm_sum)[x.sum.svd.v.genenames]


################################
# Print some genes 
################################


plot.density  <- function(xx, tt = '') {
qplot(x=xx) + cowplot::theme_cowplot() +theme(axis.title.x=element_blank(), legend.title = element_blank()) + ggtitle(tt) + scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+ scale_y_log10()
}



bal.out  <- get.excess.zeros(x_norm_sum_log, x_norm_sum_log_alra_result )

set.seed(3)
plot.density  <- function(xx, tt = '') {
#qplot(x=xx) + cowplot::theme_cowplot() +theme(axis.title.x=element_blank(), legend.title = element_blank()) + ggtitle(tt) + scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+ scale_y_log10()
qplot(x=xx[xx!=0]) + cowplot::theme_cowplot() +theme(axis.title.x=element_blank(), legend.title = element_blank()) + ggtitle(tt) + scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) + annotation_compass(sprintf('%d zeros', sum(xx==0)), 'NE') 
}
ggs  <- list()
for (i in 1:6) {
    var.features <- alra.sum.log.seurat@assays$RNA@var.features
    bal.out.var  <- bal.out[var.features]
    qindex  <-  quantile(bal.out.var, prob = c(0,0.25))
    #qindex  <-  quantile(bal.out.var, prob = c(0.75, 1))
    g1  <- var.features[bal.out.var > qindex[1] & bal.out.var < qindex[2]] %>% sample(1)
    g_raw  <- x_t[,g1]
    g_lra  <- x.lra.sum.log[g1,]
    g_alra  <- x.alra.sum.log[g1,]
    gg1 <- plot.density(g_raw) 
    gg2 <- plot.density(g_lra, sprintf('%s (%.2f)', g1, bal.out.var[g1]))
    gg3 <- plot.density(g_alra)
    gg  <- gg1 + gg2 + gg3
    ggs[[i]]  <- gg
    #print(bal.out.var[g1])
    #ggsave(gg, file=sprintf('figs/gene_histogram_%s.pdf', g1), width=8, height = 3)
}
ggg  <- patchwork::wrap_plots(ggs, ncol=1)
#ggsave(ggg, file=sprintf('figs/gene_histograms_high.pdf', g1), width=8, height = 12)
ggsave(ggg, file=sprintf('figs/gene_histograms_low.pdf', g1), width=8, height = 12)


g1  <-  'Foxd1'
bal.out[g1]

summary(g_lra)
x0 <- mt0[,g1]
x <- alra.result[[1]][,g1]
thresh  <- abs(quantile(x, 0.001))
#(sum(x0 == 0 & x > 0 & x < thresh)- sum(x0==0 & x < 0 ) ) / length(x)
(sum(x0 == 0 & x > 0 & x < thresh)- sum(x0==0 & x < 0 ) ) / sum(x0==0)
sum(x0==0 & x!=0)
sum(x0 == 0 & x > 0 & x < thresh)
sum(x0==0 & x < 0 )
length(x)







sum(g_alra == 0)

coli  <- which('Rnase6' == colnames(mt0))
bal[coli]
bal['Rnase6']
pos.vals[coli]
neg.vals[coli]


mins[coli]
g1  <- 'Foxd1'



myggplot <- ggplot(myggdata) + geom_density(aes(exp, ..density.., fill = label), alpha = 0.5) + theme_cowplot() + theme(legend.title = element_blank()) + xlab(xlab) + ylab(ylab) + ggtitle(title)


rownames(x)[9052]


################################
# Compute correlations 
################################
str(alra.sum.log.seurat@assays$RNA@data)
nonzero_idx  <- x_norm_sum_log != 0
cor(x_norm_sum_log[nonzero_idx], t(x.alra.sum.log)[nonzero_idx], method = "pearson")





################################
# What's the difference between cells super close to eachother? 
################################
#library(ranger)
#
#alra.sum.log.seurat <- RunPCA(alra.sum.log.seurat, features = VariableFeatures(alra.sum.log.seurat), verbose = F, npcs = 70)
#alra.sum.log.seurat <- RunTSNE(
#  alra.sum.log.seurat, dims = 1:70, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
#  fast_tsne_path="/home/jz437/git/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
#)
#alra.tsne <- data.frame(alra.sum.log.seurat@reductions$tsne@cell.embeddings, ident2 = ident2)
#head(alra.tsne)
#
#poi  <- 333
#nns <- FNN::knnx.index(data = alra.tsne[,1:2], query = alra.tsne[poi,1:2])
#alra.tsne.nns  <- alra.tsne[unlist(nns), 1:2]
#g <- ggplot(alra.tsne, aes(x=tSNE_1, y = tSNE_2, color = ident2)) + geom_point(shape=16, size =0.01) + theme_tsne +  geom_point(data = alra.tsne.nns, mapping=aes(x=tSNE_1, y = tSNE_2), color = 'grey', size = 0.02)
##geom_point(aes(x=tSNE_1[poi], y = tSNE_2[poi]), size = 1, color = 'black') +
#ggsave(g, file = 'figs/temp.pdf')
#
#
#dim(x_t)
#fofo <- x_t[nns[1:2],]
#fofo  <- fofo[,colSums(fofo !=0) > 1]
#dim(fofo)
#fofo %>%t
#fofo[,1:20]



################################
# BACKUP 
################################

rownames(x)[9052]
x_norm_sum_m <- colMeans(x_norm_sum)
sort(x_norm_sum_m, decreasing=T)[1:100]
x_norm_sum_m[which.max(x.sum.svd$v[,4])]
length(x_norm_sum_m)
sum(x.sum.svd$v





alra.sum.log.seurat <- RunPCA(alra.sum.log.seurat, features = VariableFeatures(alra.sum.log.seurat), verbose = F, npcs = 70)
alra.tsne <- data.frame(alra.sum.log.seurat@reductions$tsne@cell.embeddings, ident2 = ident2)
library(ranger)


alra.sum.log.seurat <- RunTSNE(
  alra.sum.log.seurat, dims = 1:70, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/home/jz437/git/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)



errs.alra <- df$Y != rfr.alra$predictions
subidx <- sample.int(nrow(alra.tsne), 2E4)
library(ggplot2)
library(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my.cols <- gg_color_hue(length(unique(ident2)))
p2 <- ggplot(alra.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2),
                                              shape=16, size =0.01)  +
  theme_tsne +
  labs(title = "ALRA")+
  scale_colour_manual(name = "Cell type", values = my.cols)  +
  annotation_compass(sprintf("Error: %.1f%%", rfr.alra$prediction.error*100), "SE")

ggsave(p2, filename = "figs/Hrvatin-tsne-alra-log-svd-nohkgenes.pdf", width = 1.45, height = 1.45)












df <- data.frame(t(x.alra.sum.log.sub[alra.sum.log.seurat@var.genes,]), Y=as.factor(ident2))
rfr.alra.sum.log <- ranger(Y ~ ., df)
errs.alra.sum.log <- df$Y != rfr.alra.sum.log$predictions










diffs <- k_choice$d[1:(length(k_choice$d)-1)] - k_choice$d[2:length(k_choice$d)] 
df <- data.frame(x=1:99,y=diffs)[10:99,]                                                 
g<-ggplot(df,aes(x=x,y=y),) + geom_point(size=2)  + geom_line()+ geom_vline(xintercept=k_choice$k)   + theme(axis.title.y=element_blank(), axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10))                                                       
g                                                                                         

ggsave('figs/spectrum_hrvatin.pdf', width=4,height=2, g)                           
print(sprintf('Chose k=%d',k_choice$k)) 

# read in dca result
fn <- "data/dca_hrvatin/mean_norm.tsv"
if (!file.exists(fn)) {
    system.time(system(sprintf('dca data/GSE102827_filtered_raw.csv data/dca_hrvatin')))
}
x.dca <- read.table(fn, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)


# magic
fn ='data/hrvatin_magic.RDS';
if ( !file.exists(fn)) {
        cat(sprintf("Generating %s\n", fn))
	system.time(	x.magic <- magic(x_norm))
        x.magic <- t(x.magic$result)
	saveRDS(x.magic,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        x.magic <- readRDS(fn)
}


# cell types
dca <- CreateSeuratObject(x.dca.sub, project = "dca")
magic <- CreateSeuratObject(x.magic.sub, project = "magic")

obs <- NormalizeData(obs)
dca@data <- log(dca@data + 1)

obs <- FindVariableGenes(obs, do.plot = FALSE)
alra.seurat <- FindVariableGenes(alra.seurat, do.plot = FALSE)
dca <- FindVariableGenes(dca, do.plot = FALSE)
magic <- FindVariableGenes(magic, do.plot = FALSE)

obs <- ScaleData(obs)
alra.seurat <- ScaleData(alra.seurat)
dca <- ScaleData(dca)
magic <- ScaleData(magic)

obs <- RunPCA(obs, pc.genes = obs@var.genes, do.print = FALSE, pcs.compute = 70)
fofo
alra.seurat <- RunPCA(alra.seurat, pc.genes = alra.seurat@var.genes, do.print = F, pcs.compute = 70)
dca <- RunPCA(dca, pc.genes = dca@var.genes, do.print = F, pcs.compute = 70)
magic <- RunPCA(magic, pc.genes = magic@var.genes, do.print = F, pcs.compute = 70)


ident3 <- factor(ident2, levels = levels(factor(ident2))[c(3, 4, 5, 6, 7, 8, 9, 1, 2, 10, 12, 13, 11)])
dim.obs <- k_choice$k
dim.dca <- k_choice$k


obs <- RunTSNE(
  obs, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/data/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)


dca <- RunTSNE(
  dca, dims.use = 1:dim.dca, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/data/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
magic <- RunTSNE(
  magic, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/data/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
dca.tsne <- data.frame(dca@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
magic.tsne <- data.frame(magic@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)


print("Done with t-SNE, beginning random forests")
library(ranger)
df <- data.frame(t(obs@scale.data[obs@var.genes,]), Y=as.factor(ident2))
rfr <- ranger(Y ~ ., df)
errs <- df$Y != rfr$predictions
df <- data.frame(t(x.alra.sub[alra.seurat@var.genes,]), Y=as.factor(ident2))
rfr.alra <- ranger(Y ~ ., df)
errs.alra <- df$Y != rfr.alra$predictions
df <- data.frame(t(x.dca.sub[dca@var.genes,]), Y=as.factor(ident2))
rfr.dca <- ranger(Y ~ ., df)
errs.dca <- df$Y != rfr.dca$predictions
df <- data.frame(t(x.magic.sub[magic@var.genes,]), Y=as.factor(ident2))
rfr.magic <- ranger(Y ~ ., df)
errs.magic <- df$Y != rfr.magic$predictions
rfr
rfr.dca
rfr.magic



subidx <- sample.int(nrow(obs.tsne), 2E4)
library(ggplot2)
library(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my.cols <- gg_color_hue(length(unique(ident2)))
p3 <- ggplot(obs.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
                                    shape=16, size =0.01) + guides(colour = FALSE) + 
                 theme_tsne + 
  labs(title = "Observed") +
  annotation_compass(sprintf("Error: %.1f%%", rfr$prediction.error*100), "SE")
p2 <- ggplot(alra.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
                                     shape=16, size =0.01)  + 
                 theme_tsne + 
                 labs(title = "ALRA")+
  scale_colour_manual(name = "Cell type", values = my.cols)  + 
  annotation_compass(sprintf("Error: %.1f%%", rfr.alra$prediction.error*100), "SE")
p4 <- ggplot(dca.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
                                    shape=16, size =0.01)  + 
                 theme_tsne + 
                 labs(title = "DCA")+
  scale_colour_manual(name = "Cell type", values = my.cols)  + 
  annotation_compass(sprintf("Error: %.1f%%", rfr.dca$prediction.error*100), "SE")
p5 <- ggplot(magic.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
                                    shape=16, size =0.01)  + 
                 theme_tsne + 
                 labs(title = "MAGIC")+
  scale_colour_manual(name = "Cell type", values = my.cols)  + 
  annotation_compass(sprintf("Error: %.1f%%", rfr.magic$prediction.error*100), "SE")
g1<-plot_grid(p3,p2,  p5 ,p4,nrow=1)
ggsave(g1,filename="figs/Hrvatin-tsne.pdf",width=5.8,height=1.45)

g1





##=============================================##
## Test SVD

set.seed(3)
system.time(k_choice <- choose_k(x_norm))
print(sprintf("k=%d was chosen", k_choice$k))
system.time(result.completed <- alra(x_norm,k_choice$k))

x.svd <- t(result.completed$A_norm_rank_k)
colnames(x.svd) <- colnames(x)
x.svd.sub <- x.svd[,!is.na(ident)]
svd.seurat <- CreateSeuratObject(x.svd.sub, project = "svd")
svd.seurat <- FindVariableGenes(svd.seurat, do.plot = FALSE)
svd.seurat <- ScaleData(svd.seurat)
svd.seurat <- RunPCA(svd.seurat, pc.genes = svd.seurat@var.genes, do.print = F, pcs.compute = 70)
dim.svd <- k_choice$k
svd.seurat <- RunTSNE(
  svd.seurat, dims.use = 1:dim.svd, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/data/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
svd.tsne <- data.frame(svd.seurat@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
df <- data.frame(t(x.svd.sub[svd.seurat@var.genes,]), Y=as.factor(ident2))
rfr.svd <- ranger(Y ~ ., df)
errs.svd <- df$Y != rfr.svd$predictions

subidx <- sample.int(nrow(svd.tsne), 2E4)
p1 <- ggplot(svd.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), 
                                             shape=16, size =0.01)  + 
  theme_tsne + 
  labs(title = "SVD")+
  scale_colour_manual(name = "Cell type", values = my.cols)  + 
  annotation_compass(sprintf("Error: %.1f%%", rfr.svd$prediction.error*100), "SE")
ggsave(p1, filename = "figs/Hrvatin-tsne-svd.pdf", width = 1.45, height = 1.45)







##=============================================##
## Test k

prediction.errors <- rep(-1, 10)
alra.test.seurats <- vector("list", 10)
#k_s <- seq(45,55,2)
k_s <- c(45,55)
#for (k_i in 2:length(k_s)) {
for (k_i in 1:2) {
    k_ <- k_s[k_i]
    result.completed <- alra(x_norm,k_) 
    x.alra.test.k <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
    colnames(x.alra.test.k) <- colnames(x)
    x.alra.test.k.sub <- x.alra.test.k[,!is.na(ident)]
    alra.test.k <- CreateSeuratObject(x.alra.test.k.sub, project = "alra.test.k")
    alra.test.k <- FindVariableGenes(alra.test.k, do.plot = FALSE)
    alra.test.k <- ScaleData(alra.test.k)
    alra.test.k<- RunPCA(alra.test.k , pc.genes = alra.test.k@var.genes, do.print = F, pcs.compute = 70)
    alra.test.k <- RunTSNE( alra.test.k, dims.use = 1:k_, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
                        fast_tsne_path="~/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12) # Seurat 2 uses FIt-SNE Version 1
    alra.test.k.tsne <- data.frame(alra.test.k@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
    alra.test.seurats[[k_i]] <- alra.test.k.tsne
    df <- data.frame(t(x.alra.test.k.sub[alra.test.k@var.genes,]), Y=as.factor(ident2))
    rfr.alra.test.k <- ranger(Y ~ ., df)
    (prediction.errors[k_i] <- rfr.alra.test.k$prediction.error)
    str(alra.test.k.tsne)
    save(list = c("prediction.errors", "alra.test.seurats"), file = "stability_backups.RData")
}
prediction.errors
k_s


subidx <- sample.int(nrow(alra.test.k.tsne), 1E4)
alra.test.k.tsne <- data.frame(alra.test.seurats[[1]], ident2 = ident2, ident3 = ident3)

p1 <- ggplot(alra.test.k.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), shape=16, size =0.01) + guides(colour = FALSE) + 
                 theme_tsne + 
  labs(title = "k=45") +
  annotation_compass(sprintf("Error: %.1f%%", prediction.errors[1]*100), "SE")

alra.test.k.tsne <- data.frame(alra.test.seurats[[6]], ident2 = ident2, ident3 = ident3)
p2 <- ggplot(alra.test.k.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), shape=16, size =0.01) + guides(colour = FALSE) + 
                 theme_tsne + 
  labs(title = "k=55") +
  annotation_compass(sprintf("Error: %.1f%%", prediction.errors[2]*100), "SE")


g1<-plot_grid(p1,p2,nrow=1)
g1

ggsave(g1,filename="figs/Hrvatin-test-k-tsne.pdf",width=3.48,height=1.74)





##=============================================##
## Test quantile

prediction.errors <- rep(-1, 3)
alra.test.seurats <- vector("list", 3)
quantile.probs <- c(0, 0.001, 0.01)
for (q_i in 1:3) {
  q.p <- quantile.probs[q_i]
  set.seed(3)
  result.completed <- alra(x_norm, k = k_choice$k, quantile.prob = q.p)
  x.alra.test.q <- as.matrix(t(result.completed$A_norm_rank_k_cor_sc))
  colnames(x.alra.test.q) <- colnames(x)
  x.alra.test.q.sub <- x.alra.test.q[,!is.na(ident)]
  alra.test.q <- CreateSeuratObject(x.alra.test.q.sub, project = "alra.test.q")
  alra.test.q <- FindVariableGenes(alra.test.q, do.plot = FALSE)
  alra.test.q <- ScaleData(alra.test.q)
  alra.test.q<- RunPCA(alra.test.q , pc.genes = alra.test.q@var.genes, do.print = F, pcs.compute = 70)
  alra.test.q <- RunTSNE( alra.test.q, dims.use = 1:k_choice$k, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
                          fast_tsne_path="/home/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12) # Seurat 2 uses FIt-SNE Version 1
  alra.test.q.tsne <- data.frame(alra.test.q@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
  alra.test.seurats[[q_i]] <- alra.test.q.tsne
  df <- data.frame(t(x.alra.test.q.sub[alra.test.q@var.genes,]), Y=as.factor(ident2))
  rfr.alra.test.q <- ranger(Y ~ ., df)
  (prediction.errors[q_i] <- rfr.alra.test.q$prediction.error)
  str(alra.test.q.tsne)
}
prediction.errors


subidx <- sample.int(nrow(alra.test.q.tsne), 1E4)
alra.test.q.tsne <- data.frame(alra.test.seurats[[1]], ident2 = ident2, ident3 = ident3)
p1 <- ggplot(alra.test.q.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), shape=16, size =0.01) + guides(colour = FALSE) + 
  theme_tsne + 
  labs(title = "quantile.prob=0") +
  annotation_compass(sprintf("Error: %.1f%%", prediction.errors[1]*100), "SE")

alra.test.q.tsne <- data.frame(alra.test.seurats[[2]], ident2 = ident2, ident3 = ident3)
p2 <- ggplot(alra.test.q.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), shape=16, size =0.01) + guides(colour = FALSE) + 
  theme_tsne + 
  labs(title = "quantile.prob=0.001") +
  annotation_compass(sprintf("Error: %.1f%%", prediction.errors[2]*100), "SE")

alra.test.q.tsne <- data.frame(alra.test.seurats[[3]], ident2 = ident2, ident3 = ident3)
p3 <- ggplot(alra.test.q.tsne[subidx,]) + geom_point(aes(x = tSNE_1, y = tSNE_2, colour = ident2), shape=16, size =0.01) + guides(colour = FALSE) + 
  theme_tsne + 
  labs(title = "quantile.prob=0.01") +
  annotation_compass(sprintf("Error: %.1f%%", prediction.errors[2]*100), "SE")


g1<-plot_grid(p1,p2,p3,nrow=1)

ggsave(g1,filename="figs/Hrvatin-test-quantile-tsne.pdf",width=5.22,height=1.74)



