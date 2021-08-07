##=======================================================##
## CITE-Seq CBMC
## ALRA
##=======================================================##


library(FNN)
library(Matrix)
library(Rmagic)
library(SAVER)
library(scImpute)
library(ggplot2)
# load Seurat v3.0
devtools::load_all("../seurat/")

source('alra_ImmGen_funs.R')

##########################
# Read in the transcriptomic and ADT data
##########################

# read in data
data_exp_rna <- read.table(
  "data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F)


# filter cells
human_ratio_sum_per_cell <- colSums(data_exp_rna[grep("HUMAN_", rownames(data_exp_rna), fixed = T),]) / colSums(data_exp_rna)
#sum(human_ratio_sum_per_cell >= 0.9)
human_cell_idx <- which(human_ratio_sum_per_cell >= 0.9)
data_exp_rna <- data_exp_rna[,human_cell_idx]


# filter mouse genes
data_exp_rna <- data_exp_rna[grep("HUMAN_", rownames(data_exp_rna), fixed = T),]
dim(data_exp_rna)


# change format
data_exp_rna <- Matrix(as.matrix(data_exp_rna), sparse = T)
mean(data_exp_rna != 0)


# get ADT data
data_exp_adt <- read.table(
  "data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
data_exp_adt <- data_exp_adt[,human_cell_idx]

mean(data_exp_adt != 0)


# normalized ADT
data_norm_adt <- read.table(
  "./data/GSE100866_CBMC_8K_13AB_10X-ADT_clr-transformed.csv.gz", 
  sep = ",", header = T, row.names = 1, stringsAsFactors = F
)
data_norm_adt <- as.matrix(data_norm_adt[,human_cell_idx])


# Seurat object
data_S <- CreateSeuratObject(counts = data_exp_rna)
data_S[["ADT"]] <- CreateAssayObject(counts = as.matrix(data_exp_adt))
data_S[["ADT"]] <- CreateAssayObject(data = as.matrix(data_norm_adt))


# normalization
data_S <- NormalizeData(data_S)
data_S@assays$ADT@data <- as.matrix(data_norm_adt)


# match RNA and ADT markers
rownames(GetAssayData(data_S, assay = "ADT"))
rownames(GetAssayData(data_S, assay = "ADT")) %in% rownames(data_S)

ADT_genes <- c( "CD3"="CD3D", "CD4"="CD4", "CD8"="CD8A", "CD45RA"="PTPRC",
               "CD56"="NCAM1", "CD16"="FCGR3A", "CD10"="MME", "CD11c"="ITGAX",
               "CD14"="CD14", "CD19"="CD19", "CD34"="CD34", "CCR5"="CCR5",
               "CCR7"="CCR7")
ADT_genes <- paste("HUMAN-", ADT_genes, sep = "")
ADT_genes %in% rownames(data_S)
names(ADT_genes) <- rownames(GetAssayData(data_S, assay = "ADT"))


# get RNA data with ADT markers
data_exp_rna_adtMarkers <- GetAssayData(data_S, assay = "RNA")[ADT_genes,]

##########################
# Impute, using various methods
##########################

data_exps <- vector("list", length=6)
names(data_exps ) <- c("original", "alra", "magic", "saver","dca","scimpute")
data_exps[['original']] <-as.matrix(GetAssayData(data_S, assay = "RNA")[ADT_genes,])

# ALRA
set.seed(3)
data_S <- RunALRA(data_S)
data_exps[['alra']] <- as.matrix(GetAssayData(data_S, assay = "alra")[ADT_genes,])

# MAGIC
fn ='data/cite_cbmc_magic.RData';
if ( !file.exists(fn)) {
        cat(sprintf("Generating %s\n", fn))
        magic_result <- magic(t(GetAssayData(data_S, assay = "RNA", slot = "data")))
	saveRDS(magic_result,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        magic_result <- readRDS(fn)
}
data_exps[['magic']] <- t(magic_result$result[,ADT_genes])



# SAVER
fn ='data/cite_cbmc_saver';
if ( !file.exists(fn)) {
        cat(sprintf("Generating %s\n", fn))
        saver_result <- saver(
          GetAssayData(data_S, assay = "RNA", slot = "counts"),
          #pred.genes = which(rownames(data_S) %in% c(ADT_genes, random_genes)), 
          pred.genes = which(rownames(data_S) %in% c(ADT_genes)), 
          pred.genes.only = T
        )
	saveRDS(saver_result,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        saver_result <- readRDS(fn)
}
data_exps[['saver']]<-saver_result$estimate

# SAVER's performance is better using the estimate (i.e. without sampling)
#data_exp_adtMarkers_saver <- sample.saver(saver_result, seed = 3)[ADT_genes,]
#data_exp_adtMarkers_saver <- log(data_exp_adtMarkers_saver + 1)
#data_exps[['saver']]<- data_exp_adtMarkers_saver

#scImpute
fn ='data/cite_cbmc_scimpute.rds';
if ( !file.exists(fn)) {
        cat(sprintf("Generating %s\n", fn))
        write.table( data_exp_rna, "data/GSE100866_CBMC_rna.csv", 
                        sep = ",", quote = F, row.names = T, col.names = T)
        scimpute( count_path = "data/GSE100866_CBMC_rna.csv",
                                        out_dir = "data/scimpute_cbmc/", 
                                        Kcluster = 12)
        scimpute_norm <- as.matrix(read.table( "data/scimpute_cbmc/scimpute_count.csv", 
                                              sep = ",", header = T, row.names = 1,
                                              stringsAsFactors = F))
        rownames(scimpute_norm) <- gsub("_","-", rownames(scimpute_norm))
        scimpute_norm <- log(normbycol(scimpute_norm, 10000) + 1)
        scimpute_norm <- Matrix(scimpute_norm, sparse = T)
        scimpute_result <- scimpute_norm[ADT_genes,]
	saveRDS(scimpute_result,fn)
}else{
        print(sprintf("Loading %s\n", fn));
        scimpute_result <- readRDS(fn)
}
data_exps[['scimpute']] <- scimpute_result



fno =  "data/cbmc_dca_out/mean_norm.tsv"
if ( !file.exists(fno)) {
    print(sprintf("Generating %s", fno))
    system('dca data/GSE100866_CBMC_rna.csv data/cbmc_dca_out')
}
dca_norm <- as.matrix(read.table( fno, sep = "\t", header = T, row.names = 1,
                                 stringsAsFactors = F))
rownames(dca_norm) <- gsub("_","-", rownames(dca_norm))
dca_norm <- log(dca_norm + 1)
data_exps[['dca']] <- dca_norm[ADT_genes,]
data_exp_adtMarkers_dca <- dca_norm[ADT_genes,]


all(colnames(GetAssayData(data_S, assay = "ADT"))==colnames( GetAssayData(data_S, assay = "RNA")))

##########################
# Check nearest neighbors
##########################
for (K in c(50,100,1000) ){
    nn.adt <- get.knn(t(GetAssayData(data_S, assay = "ADT")), K)
    nn.consistency <- rep(-1, length(data_exps))
    names(nn.consistency) <- names(data_exps)
    for (i in 1:length(data_exps) ) {
        nn.rna <- get.knn(t(data_exps[[i]]), K)
        nn.consistency[i] = mean(sapply(1:nrow(nn.adt$nn.index), function(i) length( intersect(nn.rna$nn.index[i,], nn.adt$nn.index[i,] ) ) ))
    }
    print(K)
    print(round(nn.consistency,2))
}






##########################
# Compute pairwise distances
##########################
plotCorrelationPoint <- function( x, y, annotate = T, normalize = T, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1)){
    if(normalize){
        x <- x/max(x)
        y <- y/max(y)
    }
    myggdata <- data.frame(x = x, y = y)
    myggplot <- ggplot(myggdata) + theme_cowplot() + 
        geom_point(aes(x, y), size = 0.01, alpha = 0.05) + 
        xlab(xlab) + ylab(ylab)
    if(annotate){
        mycor <- cor(x, y, use = "complete.obs", method = "spearman")
        mycortext <- sprintf("r[s]==%.2f", mycor)
        myggplot <- myggplot + annotate(x = xlim[2]*0.9, y = 0.09, label = mycortext, geom="text",parse=T) + xlim(xlim[1],xlim[2]) + ylim(ylim[1],ylim[2]) + geom_abline()
    }
    return(myggplot)
}

cell_cell_dist_mat_ADT <- as.matrix(dist(t(GetAssayData(data_S,assay = "ADT"))))
cell_cell_dist_ADT <- cell_cell_dist_mat_ADT[upper.tri(cell_cell_dist_mat_ADT, diag = F)]

# subsample points
set.seed(3)
subsample_idx <- sample( x = 1:length(cell_cell_dist_ADT), size = 50000, replace = F)
gs <- vector("list", length=6)
names(gs ) <- c("Original", "ALRA", "MAGIC", "SAVER","DCA","scImpute")
for (i in 1:length(data_exps) ) {
    dist_mat <- as.matrix(dist(t(data_exps[[i]])))
    dist_mat <-  dist_mat [upper.tri(dist_mat, diag = F)]
    gs[[i]] <- plotCorrelationPoint(
      x = cell_cell_dist_ADT[subsample_idx], y = dist_mat[subsample_idx],
      xlab = "ADT", ylab =names(gs)[i], xlim = c(0,1), ylim = c(0,1))
}

g.dist.pt <- gridExtra::arrangeGrob(grobs = gs, nrow = 2)
plot(g.dist.pt)

ggsave( filename = "figs/citeseq_cbmc_cell_dist_comparison_subsample.png", 
           plot = g.dist.pt, 
           width = 11.25, height = 7.5, units = "in")


