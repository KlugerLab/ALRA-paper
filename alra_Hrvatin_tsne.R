## Hrvatin t-sne analysis with DCA results
# SAVER-data folder can be found here: https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0 
library(ranger)
library(gridExtra)
library(ggplot2) 
library(Seurat) #Seurat Version 2
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


x_t <- t(x)
totalUMIPerCell <- rowSums(x_t);
x_norm <- sweep(x_t, 1, totalUMIPerCell, '/');
x_norm <- x_norm * 10E3
x_norm <- log(x_norm +1);

# Compute ALRA
if (!file.exists('data/hrvatin_full_alra.rds')) {
  print("Computing ALRA")
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
ct.dat <- read.csv("data/hrvatin_celltypes.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
ct.dat2$celltype[!is.na(ct.dat2$subtype)] <- 
  ct.dat2$subtype[!is.na(ct.dat2$subtype)]
table(ct.dat2$celltype)
ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]
ident <- ct.dat3[colnames(x), 4]
x.sub <- x[, !is.na(ident)]
x.dca.sub <- x.dca[, !is.na(ident)]
x.magic.sub <- x.magic[,!is.na(ident)]
x.alra.sub <- x.alra[,!is.na(ident)]
ident2 <- ct.dat3[colnames(x.sub), 4]
ident2 <- factor(ident2, levels = unique(ident2)[order(tolower(unique(ident2)))])



## Compute 2D embedding
print("Beginning the Seurat analysis")
# Seurat
obs <- CreateSeuratObject(x.sub, project = "obs")
alra.seurat <- CreateSeuratObject(x.alra.sub, project = "alra")
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
dim.alra <- k_choice$k
dim.dca <- k_choice$k


obs <- RunTSNE(
  obs, dims.use = 1:dim.obs, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
  fast_tsne_path="/data/george/Research_Local/FIt-SNE_experimental/FIt-SNE_v1/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
alra.seurat <- RunTSNE(
  alra.seurat, dims.use = 1:dim.alra, check_duplicates = FALSE, do.fast = TRUE,seed.use=3, tsne.method="FIt-SNE", 
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
alra.tsne <- data.frame(alra.seurat@dr$tsne@cell.embeddings, ident2 = ident2, ident3 = ident3)
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



