### Functions used



# Function to get only legend for a ggplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# Function to normalize count data
normbycol <- function(datain, newsum){
  dat_col <- matrix(rep(colSums(datain),nrow(datain)),byrow=T,nrow = nrow(datain))
  dataout <- (datain/dat_col) * newsum
}



# Function to evaluate zero preservation
zeroEval <- function(true.data, impute.data, true.label, impute.label){
  n.label <- length(true.label)
  
  # vector for result
  zero.ratio <- rep(0,n.label)
  names(zero.ratio) <- true.label
  non.zero.ratio <- zero.ratio
  
  for(ii in 1:n.label){
    # get data for the current label
    label.tmp <- true.label[ii]
    true.tmp <- as.numeric(true.data[,true.label == label.tmp])
    impute.tmp <- impute.data[,impute.label == label.tmp]
    
    # get genes with/without 0 count in true data
    zero.gene.tmp <- which(true.tmp <= 0)
    non.zero.gene.tmp <- which(true.tmp > 0)
    
    # get ratio of 0 in impute data of zero genes
    zero.ratio.tmp <- mean(impute.tmp[zero.gene.tmp,] <= 0)
    zero.ratio[label.tmp] <- zero.ratio.tmp
    
    # get ratio of non-0 in impute data of nonzero genes
    non.zero.ratio.tmp <- mean(impute.tmp[non.zero.gene.tmp,] > 0)
    non.zero.ratio[label.tmp] <- non.zero.ratio.tmp
  }
  
  return(cbind(zero.ratio, non.zero.ratio))
}



# Function to evaluate classification based on correlation
corEval <- function(true.data, impute.data, true.label, impute.label){
  input.n.cell <- ncol(impute.data)
  n.label <- length(true.label)
  
  # confusion matrix to record max correlation result
  cor.conf.mat <- matrix(0,n.label,n.label)
  rownames(cor.conf.mat) <- true.label
  colnames(cor.conf.mat) <- true.label
  
  # go over each cell and find the bulk profile with maximum correlation
  for(ii in 1:input.n.cell){
    all.cor.tmp <- apply(true.data, MARGIN = 2, cor, y = impute.data[,ii])
    max.label.tmp <- true.label[which.max(all.cor.tmp)]
    cor.conf.mat[impute.label[ii],max.label.tmp] <- cor.conf.mat[impute.label[ii],max.label.tmp] + 1
  }
  
  return(cor.conf.mat)
}



# Function to calculate correlation with true and other bulk profile
corPerCell <- function(true.data, impute.data, true.label, impute.label){
  input.n.cell <- ncol(impute.data)
  n.label <- length(true.label)
  
  cor.cell.res <- matrix(0, nrow = input.n.cell, ncol = 2)
  
  for(ii in 1:input.n.cell){
    # get correlation with true bulk profile
    cor.bulk.true <- cor(true.data[,true.label == impute.label[ii]], impute.data[,ii])
    
    # get max correlation with other bulk profiles
    cor.bulk.other <- max(apply(true.data[,true.label != impute.label[ii]], 2, cor, y = impute.data[,ii]))
    
    # record result
    cor.cell.res[ii,1] <- cor.bulk.true
    cor.cell.res[ii,2] <- cor.bulk.other
  }
  
  return(cor.cell.res)
}



# Functions to plot correlation results
plotCellCor <- function(cor.cell.res, input.cell.label = NULL, size = 1.5, 
                        title = "", legend = F, xlab = "", ylab = ""){
  colnames(cor.cell.res) <- c("x","y")
  
  # myggplot <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #   geom_abline(slope = 1, intercept = 0, linetype = 2)
  myggplot <- ggplot() + theme_cowplot() + geom_abline(slope = 1, intercept = 0, linetype = 2)
  
  if(is.null(input.cell.label)){
    myggplot <- myggplot + 
      geom_jitter(data = data.frame(cor.cell.res), aes(x, y), size = size, col = "blue", alpha = 0.5) + 
      xlab(xlab) + ylab(ylab) + 
      ggtitle(title) + xlim(c(0,1)) + ylim(c(0,1))
  } else {
    myggplot <- myggplot + 
      geom_jitter(data = data.frame(x = cor.cell.res[,1], y = cor.cell.res[,2], 
                                    Group = factor(input.cell.label, levels = unique(input.cell.label)[
                                      order(unique(toupper(input.cell.label)))])), 
                  aes(x, y, col = Group), size = size, alpha = 0.5) + 
      xlab(xlab) + ylab(ylab) + 
      ggtitle(title) + xlim(c(0,1)) + ylim(c(0,1)) + 
      guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))
  }
  
  if(!legend){
    myggplot <- myggplot + theme(legend.position = "none")
  } else {
    myggplot <- myggplot + theme(legend.position = "bottom", legend.title = element_blank())
  }
  
  return(myggplot)
}

# function to plot correlation results against diff
plotCellCorDiff <- function(cor.cell.res, input.cell.label = NULL, size = 1.5, 
                            title = "", legend = F, xlab = "", ylab = "", xlim = c(0,1), ylim = c(-0.1,0.4)){
  colnames(cor.cell.res) <- c("x","y")
  
  # myggplot <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #   geom_abline(slope = 1, intercept = 0, linetype = 2)
  myggplot <- ggplot() + theme_cowplot() + geom_hline(yintercept = 0, linetype = 2)
  
  if(is.null(input.cell.label)){
    myggplot <- myggplot + 
      geom_jitter(data = data.frame(cor.cell.res), aes(x, x-y), size = size, col = "blue", alpha = 0.5) + 
      xlab(xlab) + ylab(ylab) + 
      ggtitle(title) + xlim(xlim) + ylim(ylim)
  } else {
    myggplot <- myggplot + 
      geom_jitter(data = data.frame(x = cor.cell.res[,1], y = cor.cell.res[,2], 
                                    Group = factor(input.cell.label, levels = unique(input.cell.label)[
                                      order(unique(toupper(input.cell.label)))])), 
                  aes(x, x-y, col = Group), size = size, alpha = 0.5) + 
      xlab(xlab) + ylab(ylab) + 
      ggtitle(title) + xlim(xlim) + ylim(ylim) + 
      guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))
  }
  
  if(!legend){
    myggplot <- myggplot + theme(legend.position = "none")
  } else {
    myggplot <- myggplot + theme(legend.position = "bottom", legend.title = element_blank())
  }
  
  return(myggplot)
}