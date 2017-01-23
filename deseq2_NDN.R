# loading necessary packages
library(BiocParallel)
library(DESeq2)
library("pheatmap")

print('Starting R script')

register(SnowParam(7))

# reading in count data
datafilepath <- "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts/read_count_table.txt"
data <- read.table(datafilepath, header = TRUE, sep = "\t", row.names=1)

# reading in colData
coldatafilepath <- "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts/col_data.txt"
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
  
# print('Fitting DEseq model')
# 
# # fit model with DESeq
# dds <- DESeq(dds, test="LRT", reduced = ~ Day + Run, parallel = TRUE) 
# 
# # get results for change between days
# res <- results(dds, contrast = c("Day","0","3"))
# 
# # sort genes according to their significance and show most significant
# res <- res[order(res$padj),]
# head(res)
# 
# # plot MA plot
# plotMA(res, alpha=0.01, ylim=c(-10,10))
# 
# # plot volcano plot
# volcanoData <- cbind(res$log2FoldChange, -log10(res$padj))
# volcanoData <- na.omit(volcanoData)
# colnames(volcanoData) <- c("logFC", "negLogPval")
# plot(volcanoData)
