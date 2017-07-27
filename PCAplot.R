#Creating PCA plots of RNA-sequencing data
install.packages('Rcpp', dependencies = TRUE)
library("Rcpp")
install.packages('colorspace', dependencies = TRUE)
library("colorspace")
library("DESeq2")
CountTable <- read.csv("P:\\CCTI_USERS/Michelle  Miron/Lymphoid TRM paper/Analysis/RNAseq analysis/NormCounts_allDEseqBMvsLNgenes_thresholdFC1padj0.05_bldbmlln.csv", header=T, row.names=1)
head(CountTable)
integernorm_counts <- round(CountTable)
Samples <- data.frame(row.names=colnames(integernorm_counts), condition=as.factor(c(rep("bld",3),"LN-1","LN-3","LN-2","BM-1","BM-2","BM-3")))
dds <- DESeqDataSetFromMatrix(countData = integernorm_counts, colData=Samples, design=~condition)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
library(genefilter)
library("ggplot2")
#with text
#benstyle.PCA = function(object, intgroup="condition", ntop=500, x=1, y=2, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,x], PC2=pca$x[,y], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(percentVar)
  }
  
  
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", label="name")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[x] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[y] * 100),"% variance")) +
    coord_fixed() + geom_text()+theme(panel.background = element_rect(fill = 'white'))
  
}
#without text
benstyle.PCA = function(object, intgroup="condition", ntop=500, x=1, y=2, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,x], PC2=pca$x[,y], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(percentVar)
  }
  
  
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", label="name")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[x] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[y] * 100),"% variance")) +
    coord_fixed() +theme(panel.background = element_rect(fill = 'white'))
  
}
pcaData <- benstyle.PCA(rld, x=1, y=2)
pcaData