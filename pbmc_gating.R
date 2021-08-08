source("../ALRA/alra.R")
library(ggplot2)
library(cowplot)
source('convenience.R')
library(patchwork)
library(tidyverse)
if ( !file.exists('data/pbmc68k_data.rds')) {
    download.file('https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/pbmc68k_data.rds', 'data/pbmc68k_data.rds')
}
pbmc68k_data <- readRDS('data/pbmc68k_data.rds')
A  <-  pbmc68k_data$all_data[[1]]$hg19$mat
scRNA_genes  <- pbmc68k_data$all_data[[1]]$hg19$gene_symbols

colnames(A) <- scRNA_genes
num_of_genes_in_cell <- rowSums(A>0)
keep_cells = which(num_of_genes_in_cell > 400)
length(keep_cells)
A_ <- as.matrix(A[keep_cells,])

num_of_cells_in_gene <- colSums(A_>0)
keep_genes = which(num_of_cells_in_gene > 5)
length(keep_genes)
A_ <- as.matrix(A_[,keep_genes])
A_df <- as.data.frame(A_)
dim(A_)
# [1] 68579 16333


dim(A_df)
# [1] 58674 16139

A_norm <- normalize_data(A_)
dim(A_norm)



# ALRA
set.seed(3)
k_choice <- choose_k(A_norm)
system.time(alra.result <- alra(as.matrix(A_norm),k_choice$k))
A_norm_alra <- as.data.frame(alra.result[[3]])

## temp: A_norm_magic <- A_magic$result
#fn ='data/68k_pbmc_magic.RDS';
#if ( !file.exists(fn)) {
#        library(Rmagic)
#        library(reticulate)
#        use_condaenv('magic')
#        print(sprintf("Generating %s\n", fn));
#        system.time(magic_result <- magic((A_norm), verbose =T))
#	saveRDS(magic_result,fn)
#}else{
#        print(sprintf("Loading %s\n", fn));
#        magic_result <- readRDS(fn)
#}
#A_magic <- magic_result$result





################################
# Gating 
################################

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
    pos_marks <- rowSums(A_norm_alra[, celltypes[[celltype]]$pos_marks] > 0  ) == length(celltypes[[celltype]]$pos_marks)
    if ( length(celltypes[[celltype]]$neg_marks) >0 ){
        neg_marks <- rowSums(A_norm_alra[, celltypes[[celltype]]$neg_marks] > 0  ) ==0
        identified <- pos_marks & neg_marks
    }else{
        identified <- pos_marks
    }
    cell_labels[identified] <- celltype
    print(celltype)
    print(sum(identified))
}
A_df$cell_labels <- cell_labels
colnames(A_df)  <- make.unique(colnames(A_df))

for ( li in 1:length(celltypes)) {
    cat(sprintf('%s & %s & %s\\\\\n',names(celltypes)[li], paste(celltypes[[li]]$pos_marks,collapse = ', '),  paste(celltypes[[li]]$neg_marks,collapse = ', ')))
}


################################
# Visualize the gating 
################################

# Run PCA using SVD. 
col_gene_means <- colMeans(A_norm)
A_norm_c <- sweep(A_norm, 2, col_gene_means, '-')
system.time(svdout <- rsvd(A_norm_c, k=k_choice$k))
obsPCA <- svdout$u %*% diag(svdout$d) 


# Run t-SNE
source('~/Research_Local/FIt-SNE/fast_tsne.R', chdir=T)
init <- 0.0001*(obsPCA[,1:2]/sd(obsPCA[,1]))
fitsneout <- fftRtsne( obsPCA,initialization = init , rand_seed=3)

#dframe <- dframe[!is.na(dframe$Y),] # Remove cells that don't have an assigned celltype
dframe1=data.frame(fitsneout, Y= as.factor(cell_labels))
(g1 <- ggplot(dframe1) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+theme(legend.position="bottom")  + guides(colour = guide_legend(nrow=3,override.aes = list(size=10)))  +  theme_cowplot() + theme(legend.position= "bottom", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 10, face = "plain"), plot.margin=unit(c(0.05,0.05,0.05,0.05 ), "cm"), axis.line.y.left=element_line(colour="black",size=0.2), axis.line.x.bottom=element_line(colour="black",size=0.2)) )
#+ labs(x="t-SNE 1", y="t-SNE 2")
ggsave(g1,filename="figs/gated_68k.pdf")

################################
# Backup 
################################

################################
# Chosen cell types 
################################
#
#blood_cell_gene_data <- read_tsv('data/rna_blood_cell.tsv')
#unique( blood_cell_gene_data$`Blood cell`)
#blood_cell_gene_data <- blood_cell_gene_data %>% rename (gene_name = `Gene name`, blood_cell = `Blood cell`) %>% select(gene_name,blood_cell, TPM) %>% pivot_wider(id_cols = blood_cell, names_from =gene_name, values_from =TPM, values_fn = list(TPM = mean   ))
#blood_cell_gene_data <- blood_cell_gene_data %>% mutate(blood.cell.type.chosen = case_when( blood_cell == 'classical monocyte' ~ 'classical_monocyte',
#                                                                   blood_cell  == 'myeloid DC' ~'myeloidDC',
#                                                                   blood_cell  == 'plasmacytoid DC' ~'plasmacytoidDC',
#                                                                   blood_cell == 'T-reg' ~ 'Treg',
#                                                                   blood_cell == 'memory CD4 T-cell' ~ 'memoryCD4',
#                                                                   blood_cell == 'memory CD8 T-cell' ~ 'memoryCD8',
#                                                                   blood_cell == 'memory B-cell' ~ 'memoryB',
#                                                                   blood_cell == 'naive B-cell' ~ 'naiveB',
#                                                                   blood_cell == 'NK-cell' ~ 'nkCells',
#                                                                   blood_cell == 'classical monocyte' ~ 'classical_monocyte',
#                                                                   blood_cell == 'intermediate monocyte' ~ 'intermediate_monocyte',
#                                                                   blood_cell == 'non-classical monocyte' ~ 'nonclassical_monocyte',
#                                                                   T ~ NA_character_))  %>% 
#                                                filter( !is.na(blood.cell.type.chosen))
#
#blood_cell_gene_data %>% group_by(blood.cell.type.chosen) %>% summarize( mean(CD27), mean(CD19), mean(CD14), mean(CD3G)  )
#blood_cell_gene_data %>% group_by(blood.cell.type.chosen) %>% summarize( mean(CD19),  mean(CD14), mean(CD3G), mean(CD27)  )
#with(blood_cell_gene_data, table(blood_cell, blood.cell.type.chosen , useNA="ifany"))
#
#
#
## Identify genes in common
#common_gn <- intersect(colnames(A_df) , colnames(blood_cell_gene_data))
## Identify variable genes
##compute_cv <- function(x) sd(x) / mean(x)
##cvs <- apply( A_df[,common_gn],2, FUN=compute_cv)
##common_gn_var <- common_gn[which(cvs > quantile(cvs,0.7))]
## Compute fcorrelations
##A_df_corrs <- cor( A_df[,common_gn_var] %>% t, blood_cell_gene_data[, common_gn_var] %>%t)
##A_df$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(A_df_corrs, 1, which.max)]
##with(A_df, table(cell_labels, most.correlated.cell.type, useNA="ifany"))
#common_types <- intersect( A_df$cell_labels, (blood_cell_gene_data$blood.cell.type.chosen))
#A_df_corrs <- cor( A_df[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t, method = 'spearman')
#A_df$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(A_df_corrs, 1, which.max)]
#tbl_ <- with(A_df, table(cell_labels, most.correlated.cell.type, useNA="ifany"))
#tbl_[c(common_types, 'none'),common_types]
#
###Try with ALRA
##A_norm_corrs <- cor( A_norm[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t)
##A_norm$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(A_norm_corrs, 1, which.max)]
##with(A_norm, table(cell_labels, most.correlated.cell.type, useNA="ifany"))
##Try with ALRA
##A_norm_alra_corrs <- cor( A_norm_alra[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t)
##A_norm_alra$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(A_norm_alra_corrs, 1, which.max)]
##with(A_norm_alra, table(cell_labels, most.correlated.cell.type, useNA="ifany"))
#
#
#
#
#################################
## Mean Esxpression of specific markers 
#################################
#A_df %>% group_by(cell_labels) %>% summarize_at( vars(FMO1, KCNG1, CXorf65, DSP, TAS1R3, FANK1, SEMA3G, S1PR3, ANXA3, CRIP3), function(r) mean())
#
##A_ %>% as.matrix %>% as.data.frame %>% setNames(make.names(names(.), unique = TRUE))  %>%  summarize_at( vars(FMO1, KCNG1, CXorf65, DSP, TAS1R3, FANK1, SEMA3G, S1PR3, ANXA3, CRIP3), function(r) sum(r>0), )
##A_ %>% as.matrix %>% as.data.frame %>% setNames(make.names(names(.), unique = TRUE))  %>%  summarize_at( vars(TAS1R3,  FA2H, PPP1R37,     SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) sum(r>0), )
#
#
#
#A_temp <- as.data.frame(as.matrix(A[keep_cells,])) %>% setNames(make.names(names(.), unique = TRUE)) 
#A_temp$cell_labels <- A_df$cell_labels
#A_temp[1:5,1:5]
#table(A_temp$cell_labels)
#A_temp  %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) sum(r>0))
#blood_cell_gene_data  %>%group_by(blood.cell.type.chosen) %>% select( TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1) # In the bulk RNA-seq
#table( A_df$cell_labels[A_temp$PPP1R37 >0 ], useNA="ifany") %>% t %>% t
#
#A_temp %>% group_by(cell_labels) %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) mean(r))
#A_temp %>% group_by(cell_labels) %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37,    OR4C6, SDC1, EYA2, ARL14, TMEM74,AMPD1,TUBB2B,SAA4,PROX1), function(r) mean(r))
#
#A_temp %>%setNames(make.names(names(.), unique = TRUE)) %>% group_by(A_df$most.correlated.cell.type) %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37   ), function(r) mean(r>0))
#A_norm_alra %>%setNames(make.names(names(.), unique = TRUE)) %>% group_by(A_df$cell_labels) %>% summarize_at( vars(TAS1R3,  FA2H, PPP1R37   ), function(r) mean(r>0))
#
##
#A_df %>% group_by(cell_labels) %>% summarize(n())
##A_df %>% group_by(most.correlated.cell.type) %>% summarize(n())
#
##A_df %>% group_by(cell_labels) %>% summarize_at( vars(CD3G, CD4, CD8A, PTPRC), function(r) mean(r>0))
##
##'HLA-DRA' %in% colnames(A_df)
##A_df %>% group_by(cell_labels) %>% summarize( FMO1=mean(FMO1), FMO1=mean(FMO1)  )
##A_df %>% group_by(cell_labels) %>% summarize_at( vars(FMO1, KCNG1, CXorf65, DSP, TAS1R3, FANK1, SEMA3G, S1PR3, ANXA3, CRIP3), function(r) mean(r>0))
##                                pos_marks = c('FUT4','FCGR3A'),
##                                neg_marks = c('CD19', 'CD19')
##A_df %>% group_by(cell_labels) %>% summarize_at( vars(FUT4, FCGR3A, CD19), function(r) mean(r>0))
##A_df %>% group_by(most.correlated.cell.type) %>% summarize_at( vars(FUT4, FCGR3A, CD19), function(r) mean(r>0))
##A_df %>% group_by(most.correlated.cell.type) %>% summarize_at( vars(CD3G, CD8A, CD4, PTPRC), function(r) mean(r>0))
##A_df %>% filter(most.correlated.cell.type == 'memoryCD8' & cell_labels == 'neutrophils') %>% summarize_at( vars(FUT4, FCGR3A,CD19), function(r) mean(r>0))
##A_df %>% filter(most.correlated.cell.type == 'memoryCD8' & cell_labels == 'neutrophils') %>% summarize_at( vars(CD3G, CD8A, CD4, PTPRC), function(r) mean(r>0))
##blood_cell_gene_data[,c('blood_cell','FUT4', 'FCGR3A','CD19')]
##blood_cell_gene_data[,c('blood_cell','PTPRC', 'FCGR3B', 'FCGR3A','CD19')]
##blood_cell_gene_data[,c('blood_cell','CD3G', 'CD8A', 'CD4', 'PTPRC')]
##
##table( A_df$cell_labels, useNA="ifany")
##tail(colnames(A_df))
#
#
#
#
#
#################################
## Visualize! 
#################################
#
## Run PCA using SVD. 
#col_gene_means <- colMeans(A_norm)
#A_norm_c <- sweep(A_norm, 2, col_gene_means, '-')
#system.time(svdout <- rsvd(A_norm_c, k=k_choice$k))
#obsPCA <- svdout$u %*% diag(svdout$d) 
#
#
## Run t-SNE
#source('~/Research_Local/FIt-SNE/fast_tsne.R', chdir=T)
#init <- 0.0001*(obsPCA[,1:2]/sd(obsPCA[,1]))
#fitsneout <- fftRtsne( obsPCA,initialization = init , rand_seed=3)
#
##dframe <- dframe[!is.na(dframe$Y),] # Remove cells that don't have an assigned celltype
#dframe1=data.frame(fitsneout, Y= as.factor(cell_labels))
#levels(dframe$Y)
#(g1 <- ggplot(dframe1) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+theme(legend.position="bottom") + labs(x="t-SNE 1", y="t-SNE 2") + guides(colour = guide_legend(override.aes = list(size=10))) + ggtitle('Labled by markers'))
#dframe=data.frame(fitsneout, Y= as.factor(A_df$most.correlated.cell.type ))
#dframe$Y  <- factor(dframe$Y, levels=c( levels(dframe1$Y), setdiff(levels(dframe$Y), levels(dframe1$Y))))
#(g2 <- ggplot(dframe) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+theme(legend.position="bottom") + labs(x="t-SNE 1", y="t-SNE 2") + guides(colour = guide_legend(override.aes = list(size=10))) + ggtitle('Labeled by correlations'))
#ggsave(g1+g2, file="figs/Zheng_tsne_labels.pdf", width=14)
#
#
##dframe <- dframe[!is.na(dframe$Y),] # Remove cells that don't have an assigned celltype
#dframe1=data.frame(fitsneout, Y= as.factor(cell_labels)) %>% filter( Y %in% common_types)
#dframe1=data.frame(fitsneout, Y= as.factor(cell_labels)) %>% filter( Y %in% c('memoryCD4', 'memoryCD8'))
#levels(dframe$Y)
#(g1 <- ggplot(dframe1) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+theme(legend.position="bottom") + labs(x="t-SNE 1", y="t-SNE 2") + guides(colour = guide_legend(override.aes = list(size=10))) + ggtitle('Labled by markers'))
#dframe=data.frame(fitsneout, Y= as.factor(A_df$most.correlated.cell.type )) %>% filter( Y %in% common_types)
#dframe=data.frame(fitsneout, Y= as.factor(A_df$most.correlated.cell.type )) %>% filter( Y %in% c('memoryCD4', 'memoryCD8'))
#dframe$Y  <- factor(dframe$Y, levels=c( levels(dframe1$Y), setdiff(levels(dframe$Y), levels(dframe1$Y))))
#(g2 <- ggplot(dframe) + geom_point(aes(x=X1, y=X2, color=Y), size=0.5)+theme(legend.position="bottom") + labs(x="t-SNE 1", y="t-SNE 2") + guides(colour = guide_legend(override.aes = list(size=10))) + ggtitle('Labeled by correlations'))
#ggsave(g1+g2, file="figs/Zheng_tsne_labels.pdf", width=14)
#
#
#################################
## BU 
#################################
#
#
#
#
#
#
#
#A_df$num_of_nonzeros  <- rowSums(A_df[,common_gn] != 0)
#summary(A_df$num_of_nonzeros)
#
#
#blood_cell_gene_data$num_non_zero <- rowSums(blood_cell_gene_data[,common_gn] > 0)
#blood_cell_gene_data %>% select(num_non_zero, blood.cell.type.chosen)
#
#
#
#length(common_gn)
#
#
#
#
#
#
#
#A_df_common <- A_df[A_df$cell_labels != "none", ] %>% as_tibble
#cvs <- apply( blood_cell_gene_data[,common_gn], 2, compute_cv)
#common_gn_variable <- common_gn[which(cvs > 1.2)]
#
#fofo <- cor( a_df_common[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t, method = 'spearman')
#colnames(fofo)  <-  blood_cell_gene_data$blood_cell
#toplot <- fofo[A_df_common$cell_labels == 'neutrophils',] %>% as_tibble  %>% pivot_longer( everything(),names_to="bulk_cell_type", values_to="expression")
#g <- ggplot(toplot, aes(color=bulk_cell_type, x=expression)) + geom_density()
#ggsave(sprintf('figs/temp.pdf'), g)
#
#
#
#
#
#
#fofo <- cor( A_df[,common_gn] %>% t, blood_cell_gene_data[, common_gn] %>%t, method = 'spearman')
#colnames(fofo)  <-  blood_cell_gene_data$blood_cell
#
#
#
#
#
#
#
#
#
#
#A_df_common$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(fofo, 1, which.max)]
#A_df_common %>% select( cell_labels, most.correlated.cell.type)
#
#table( A_df_common$cell_labels, useNA="ifany")
#table( A_df_common$cell_labels, A_df_common$most.correlated.cell.type, useNA="ifany")
#mean(A_df_common$cell_labels == A_df_common$most.correlated.cell.type)
#
#
#
#A_norm_alra$cell_labels <- A_df$cell_labels
#A_norm_alra_common <- A_norm_alra[A_norm_alra$cell_labels != "none", ]
#
#fofo <- cor( A_norm_alra_common[,common_gn_variable] %>% t, blood_cell_gene_data[, common_gn_variable] %>%t, method = 'spearman')
#A_norm_alra_common$most.correlated.cell.type <- blood_cell_gene_data$blood.cell.type.chosen[ apply(fofo, 1, which.max)]
#table( A_norm_alra_common$cell_labels, A_norm_alra_common$most.correlated.cell.type, useNA="ifany")
#mean(A_norm_alra_common$cell_labels == A_norm_alra_common$most.correlated.cell.type)
#
#
#compute_cv <- function(x) sd(x) / mean(x)
#
#
## Do they even match eachother?
#
#fifi <- cor(  blood_cell_gene_data[, common_gn_variable] %>%t, method = 'spearman')
#colnames(fifi)  <- blood_cell_gene_data$blood.cell.type.chosen
#rownames(fifi)  <- blood_cell_gene_data$blood.cell.type.chosen
#fifi
#
#A <- scale(blood_cell_gene_data[, common_gn], scale = F)
#svdout <- svd(A)
#g <- qplot(x=svdout$u[,1], y=svdout$u[,2]) + geom_point()
#ggsave(sprintf('figs/temp.pdf'), g)
#cor(svdout$u[,1:3] %>%t)
#svdout$u[,1:2]
#svdout$u
#svdout$d
#
#
#
#dim(fofo)
#
#
#library(corrr)
#fofo  <-  A_df_common[, common_gn] %>% correlate( blood_cell_gene_data[, common_gn] )
#
#fofo$cell_labels  <- A_df_common$cell_labels
#dim(fofo)
#
#blood_cell_gene_data[,1:5]
#
#
#
#A_df[1:10, common_gn]
#blood_cell_gene_data[, common_gn]
#A_gn_intersect <- A
#
#
#
#
#
#
#singlecell_blood[1:5,1:5]
#singlecell_blood <-  A_df %>% group_by(cell_labels) %>% summarize_all(  mean)
#A_ave <- as.data.frame(t(singlecell_blood)[-1,])
#colnames(A_ave)  <- singlecell_blood$cell_labels
#A_ave$gene_names  <- colnames(singlecell_blood)[-1]
#
#hpa_blood[1:5,1:5]
#joined <- inner_join(hpa_blood, A_ave, by = c('Gene name'='gene_names'))
#dim(hpa_blood)
#dim(A_ave)
#dim(joined)
#colnames(joined)
#joined2 <- sapply(joined, function(x) as.numeric(as.character(x)))
#joined[1:10,c(112:119)]
#joined2[1:10,c(112:119)]
#cor( joined2[,c(3:111)], joined2[,c(112:119)])
#cor( joined[,c(3:4)], joined[,c(112:119)])
#joined[,c(3:111)]
#joined[,c(112:119)]
#A_ave[16135:16137,]
#singlecell_blood[1:5,1:5]
#
#
#
#hpa_blood <- read_tsv('data/rna_blood_cell_sample_tpm_m.tsv')
#dim(hpa_blood)
#
#
#
#
#cell_labels <- rep('none', nrow(A_norm))
#for ( celltype in names(celltypes)) {
#    pos_marks <- rowSums(A_magic[, celltypes[[celltype]]$pos_marks] > 0  ) == length(celltypes[[celltype]]$pos_marks)
#    if ( length(celltypes[[celltype]]$neg_marks) >0 ){
#        neg_marks <- rowSums(A_magic[, celltypes[[celltype]]$neg_marks] > 0  ) ==0
#        identified <- pos_marks & neg_marks
#    }else{
#        identified <- pos_marks
#    }
#    cell_labels[identified] <- celltype
#    print(celltype)
#    print(sum(identified))
#}
#
#
#
#
#
#
#library(ggplot2)
#library(reshape2)
#unlist(unlist(celltypes))
#
#
#A_magic[1:5,1:5]
#
#class(A_magic)
#
#A_pos  <- melt(A_norm_alra[, celltypes$memoryCD4$pos_marks])
#A_neg  <- melt(A_norm_alra[, celltypes$memoryCD4$neg_marks])
#
#celltype
#for ( celltype in names(celltypes)) {
#    A_pos  <- melt(A_magic[, celltypes[[celltype]]$pos_marks])
#    A_neg  <- melt(A_magic[, celltypes[[celltype]]$neg_marks])
#    #A_posneg <- bind_rows(list('positive'=A_pos, 'negative'=A_neg), .id = 'type')
#    g1  <- ggplot(A_pos,aes(x=value)) +  geom_density(aes(y=..scaled..), alpha=1/2)   + facet_wrap(~variable, ncol=3)
#    g2  <- ggplot(A_neg,aes(x=value)) +  geom_density(aes(y=..scaled..), alpha=1/2)   + facet_wrap(~variable, ncol=3)
#    g <- g1 / g2
#    ggsave(sprintf('figs/magic_%s_markers.pdf', celltype), g)
#}
#
#
#
#
#
#
##g  <- ggplot(A_toplot %>% filter(variable == 'CD3G'),aes(x=value)) +  geom_density(aes(y=..scaled..), alpha=1/2)  
##p <- ggplot_build(g)
##str(p$data[[1]]$density)
##dens <- -p$data[[1]]$density
#
##fifi <- pracma::findpeaks(dens)[1,1]
