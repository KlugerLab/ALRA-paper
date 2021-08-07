# Function to overlay average expression of a certain gene in clusters
ClstAvgExp_facetted <- function(
  normData, label, labs, genes, geneGroup,methods,
  lowCol = "blue", highCol = "red", title = "", scale = 5, return.data = F
){
  actual_gene_names <- unlist(lapply( strsplit(rownames(normData),split = ".", fixed=T), function(x) x[2]))
  subData <- normData[order(actual_gene_names),]
  subLabel <- label[label %in% labs]
  
  
  avgData <- apply(subData, 1, FUN = function(x, xlabels){
    by(data = x, INDICES = xlabels, FUN = sum, na.rm = T) /
      (by(data = x != 0.00, INDICES = xlabels, FUN = sum, na.rm = T) + 0.0001)
  },xlabels = subLabel)
  avgData[avgData <0] <- 0
  avgData <- sweep(avgData,2, apply(avgData,2,max, na.rm = T), '/')
  
  avgData <- melt(avgData, varnames = c("Cluster", "Gene"), value.name = "AverageExpression")
  numData <- melt(apply(subData, 1, FUN = function(x, xlabels){
    by(data = x > 0.00, INDICES = xlabels, FUN = sum, na.rm = T) / table(xlabels) * 100
    #by(data = x, INDICES = xlabels, FUN = length)
  }, xlabels = subLabel), varnames = c("Cluster", "Gene"), value.name = "NumCells")
  
  myggdata <- cbind(avgData, numData[,3])
  # print(unique(myggdata$Cluster))
  colnames(myggdata)[4] <- "NumCells"
  
  myggdata$genegroup <- factor(unlist(lapply( strsplit(as.character(myggdata$Gene),split = ".", fixed=T), 
                                              function(x) x[2])), levels=genes)
  myggdata$method <- (unlist(lapply( strsplit(as.character(myggdata$Gene),split = ".", fixed=T), function(x) x[1])))
  myggdata$Cluster<- factor(myggdata$Cluster,levels=c(labs))
  myggdata$method<- factor(myggdata$method,levels=methods)
  #myggdata$genegroup <- factor(actual_gene_names)
  #names(geneGroup) <- genes
  #myggdata$genegroup <- geneGroup[myggdata$Gene]
  
  if(return.data){
    return(myggdata)
  }
  
  myggplot <- ggplot(data = myggdata ,aes(x = method, y = (Cluster), size = NumCells, col = AverageExpression)) + 
    geom_point() + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_grid(.~genegroup,scales="free_x") + 
    ylab("Cluster") + scale_colour_gradient(low = lowCol, high = highCol, guide_colourbar(title = "Average\nExpression\n")) + 
    scale_size(name = "% Cells > 0", range = c(0, scale)) + 
    guides(size = guide_legend(order=1), color = guide_colourbar(order = 2)) + ggtitle(title) + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
          axis.title.y=element_blank(),axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  (myggplot)
}