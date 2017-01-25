## Hey Esther, 
## Ik heb dit R script een beetje aangepast zodat het zou moeten werken op de data die ik je mee heb gegeven.
## Hetscript zoekt de data files in de folder waar het zich in bevindt, 
## dus vergeet deze niet mee in dezelfde map te steken.
## Om dit script te laten werken moet je een aantal packages installeren.
## Om DESeq2 te installeren:
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
## Om de heatmap package te installeren
install.packages("pheatmap")
## Je hoeft deze lijnen maar een keer uit te voeren en kan ze daarna verwijderen.
## Ik heb hier en daar wat extra uitleg geschreven als comments met ## voor,
## de comments met een enkele # stonden al in het script en beschrijven kort wat de onderstaande code doet.
## Veel succes met de analyse!

# loading necessary packages
library(DESeq2)
library("pheatmap")

# reading in count data
datafilepath <- "combined_lane_counts.txt"
data <- read.table(datafilepath, header = TRUE, sep = "\t", row.names=1)

# reading in colData
coldatafilepath <- "colData.txt"
colData <- read.table(coldatafilepath, header = TRUE, sep = "\t")
colData$Ind <- as.factor(colData$Ind)
colData$Day <- as.factor(colData$Day)
colData$Repeat <- as.factor(colData$Repeat)
colData$Run <- as.factor(colData$Run)

# number of patients included in the data
num_pats <- length(unique(colData$Ind))

# create DESeq data
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ Ind + Day + Run)

# plot PCA
rld <- varianceStabilizingTransformation(dds, blind = TRUE) ## deze stap kan lang duren!
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plotPCA(rld, intgroup = c("Run")) ## je kan "Ind" veranderen naar Day of Run om daar de plot van te zien

# plot heatmap
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Day, sep="-" )
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)

# filter out genes with low read counts accross samples
dds <- dds[ rowSums(counts(dds)) > 6*num_pats, ]

# fit model with DESeq
dds <- DESeq(dds, test="LRT", reduced = ~ Day + Run) ## deze stap kan lang duren!

## voor het verschil tussen twee dagen te zien moet je de eerstvolgende lijn aanpassen
## afhankelijk van welk verschil je wilt zien.
## Voor verschil dag 0 - dag 3:
## res <- results(dds, contrast = c("Day","0","3"))
## Voor verschil dag 0 - dag 7:
## res <- results(dds, contrast = c("Day","0","7"))
## Voor verschil dag 3 - dag 7:
## res <- results(dds, contrast = c("Day","3","7"))
## (de rest van de code blijft hetzelfde)

# get results for change between days
res <- results(dds, contrast = c("Day","0","3"))

# sort genes according to their significance and show most significant
res <- res[order(res$padj),]
head(res)

# plot MA plot
plotMA(res, alpha=0.01, ylim=c(-10,10))

# plot volcano plot
volcanoData <- cbind(res$log2FoldChange, -log10(res$padj))
volcanoData <- na.omit(volcanoData)
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData)
