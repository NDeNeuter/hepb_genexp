#!/usr/bin/R
library(DESeq2)

setwd("/Users/pietermeysman/Data/Immunology/GOA/HepB gen exp/reads/")

ILL=read.table('allcounts.txt', sep = '\t', header=T, row.names=1)
head(ILL)

colData = read.table('sampledesign.txt',header = T, sep = "\t", row.names=1)
colData$Run = factor(colData$Run)

dds <- DESeqDataSetFromMatrix(countData = ILL, colData = colData, design = ~ Patient + Time + Run)

dds <- dds[ rowSums(counts(dds)) > 100, ]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")

as.data.frame(colData(dds))
resultsNames(dds)

res03 <- results(dds, contrast=list('TimeEXP0','TimeEXP3'))
res03
DESeq2::plotMA(res03, ylim=c(-2,2))
res37 <- results(dds, contrast=list('TimeEXP3','TimeEXP7'))
res37
DESeq2::plotMA(res37, ylim=c(-2,2))
res07 <- results(dds, contrast=list('TimeEXP0','TimeEXP7'))
res07
DESeq2::plotMA(res07, ylim=c(-2,2))
resp12 <- results(dds, contrast=list('PatientH1','PatientH2'))
resp12


write.csv(as.data.frame(res03)[order(res03$pvalue),],
          file="RNAseq_mixed_day0vs3.csv")
write.csv(as.data.frame(res07)[order(res07$pvalue),],
          file="RNAseq_mixed_day0vs7.csv")
write.csv(as.data.frame(res37)[order(res37$pvalue),],
          file="RNAseq_mixed_day3vs7.csv")
write.csv(as.data.frame(resp12)[order(resp12$pvalue),],
          file="RNAseq_mixed_p1vs2.csv")


se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),colData=colData(dds))
plotPCA(DESeqTransform(se), intgroup=c("Time","Run"),ntop=3000,returnData=FALSE)
pcase<-plotPCA(DESeqTransform(se), intgroup=c("Time","Run","Patient"),ntop=10000,returnData=TRUE)

require('ggplot2')
ggplot(pcase,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Run))
ggplot(pcase,aes(x=PC1,y=PC2)) + geom_point(aes(colour=Time))
ggplot(pcase,aes(x=PC1,y=PC2,label=colData$Patient)) + geom_text(aes(colour=colData$Patient))
