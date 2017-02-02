# loading necessary packages
library(BiocParallel)
library(DESeq2)
library("pheatmap")
library(ggplot2)

register(SnowParam(7))

# reading in count data
datafilepath <- "read_count_table.txt"
data <- read.table(datafilepath, header = TRUE, sep = "\t", row.names=1)

# reading in colData
coldatafilepath <- "col_data.txt"
colData <- read.table(coldatafilepath, header = TRUE, sep = "\t")
colData$Ind <- as.factor(colData$Ind)
colData$Day <- as.factor(colData$Day)
colData$Repeat <- as.factor(colData$Repeat)
colData$Resp <- as.factor(colData$Resp)
colData$Run <- as.factor(colData$Run)

# number of patients included in the data
num_pats <- length(unique(colData$Ind))

# create DESeq data
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ Resp + Day + Run)

# filter out genes with low read counts accross samples
dds <- dds[ rowSums(counts(dds)) > 6*num_pats, ]

# plot PCA
rld <- varianceStabilizingTransformation(dds, blind = TRUE)
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
pdf(file = paste(getwd(),'/pcaplot.pdf', sep = ''))
plotPCA(rld, intgroup = c("Day"))
dev.off()

# plot heatmap
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Day, sep="-" )
pdf(file = paste(getwd(),'/heatmap.pdf', sep = ''))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
dev.off() 

# fit model with DESeq
dds <- DESeq(dds, test="LRT", reduced = ~ Day + Run, parallel = TRUE)

for (x in list(c("Day", 0, 3), c("Day", 0, 7), c("Day", 3, 7), c("Resp", "Non-res", "Resp"))){
  
  # get results for change between days
  res <- results(dds, contrast = c(x[1], x[2], x[3]))
  
  # sort genes according to their significance and show most significant
  res <- res[order(res$padj),]
  write.csv(as.data.frame(res), file = paste(getwd(),'/results_',x[1],'_',x[2],'.txt', sep = ''))
  head(res)
  
  # plot MA plot
  pdf(file = paste(getwd(),'/MAplot_',x[1],'_',x[2],'.pdf', sep = ''))
  plotMA(res, alpha=0.05, ylim=c(-10,10))
  dev.off()
  
  # plot volcano plot
  volcanoData <- cbind(res$log2FoldChange, -log10(res$padj))
  volcanoData <- na.omit(volcanoData)
  colnames(volcanoData) <- c("10logFC", "negLogPval")
  pdf(file = paste(getwd(),'/volcanoplot_',x[1],'_',x[2],'.pdf', sep = ''))
  plot(volcanoData)
  dev.off()
  
}