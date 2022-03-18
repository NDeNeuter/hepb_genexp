#!/usr/bin/R
library(DESeq2)

ILL<-read.table('../data/readcounts.txt', sep = '\t', header=T, row.names=1)
head(ILL)

colData <- read.table('../data/sampledesign.txt',header = T, sep = "\t", row.names=1)
colData$Response <- factor(colData$Response)

# DEGs between time points for all

dds <- DESeqDataSetFromMatrix(countData = ILL, colData = colData, design = ~ Patient + Time)

dds <- dds[ rowSums(counts(dds)) > 100, ]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")

as.data.frame(colData(dds))
resultsNames(dds)

res3 <- results(dds, contrast=list("TimeEXP3","TimeEXP0"))
DESeq2::plotMA(res3, ylim=c(-2,2))
write.csv(as.data.frame(res3)[order(res3$pvalue),],
          file="../DEG_output/All-EXP3-EXP0.csv")

res7 <- results(dds, contrast=list("TimeEXP7","TimeEXP0"))
DESeq2::plotMA(res7, ylim=c(-2,2))
write.csv(as.data.frame(res7)[order(res7$pvalue),],
          file="../DEG_output/All-EXP7-EXP0.csv")

# DEGs between responder types

dds <- DESeqDataSetFromMatrix(countData = ILL, colData = colData, design = ~ Response + Time)

dds <- dds[ rowSums(counts(dds)) > 100, ]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")

as.data.frame(colData(dds))
resultsNames(dds)

res <- results(dds, contrast=list('ResponseNR','ResponseR'), independentFiltering=TRUE)
DESeq2::plotMA(res, ylim=c(-2,2))
write.csv(as.data.frame(res)[order(res$pvalue),],
          file="../DEG_output/R-vs-NR.csv")

res3 <- results(dds, contrast=list("ResponseR.TimeEXP0","ResponseNR.TimeEXP0"))
DESeq2::plotMA(res3, ylim=c(-5,5))
write.csv(as.data.frame(res3)[order(res3$pvalue),],
          file="../DEG_output/R-vs-NR-EXP0.csv")



se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),colData=colData(dds))
plotPCA(DESeqTransform(se), intgroup=c("Response"),ntop=3000,returnData=FALSE)
plotPCA(DESeqTransform(se), intgroup=c("Patient"),ntop=3000,returnData=FALSE)
pcase<-plotPCA(DESeqTransform(se), intgroup=c("Time","Patient"),ntop=10000,returnData=TRUE)

require('ggplot2')
ggplot(pcase,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Time))

ggplot(pcase,aes(x=PC1,y=PC2,label=colData$Patient)) + geom_text(aes(colour=colData$Time)) +
  xlab("PC1 (9% Variance)") + ylab("PC2 (2% Variance)") + 
  scale_colour_discrete(name="Time Point",
                      breaks=c("EXP0", "EXP3", "EXP7"),
                      labels=c("Day 0", "Day 3", "Day 7"))
