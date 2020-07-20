# cleaning
rm(list = ls())

library(readxl)
library(dplyr)
library(DESeq2)

# loading data
setwd("~/Desktop/HT projects/VPS37_RNASeq")

# load data for DESeq2 analysis
cts <- as.matrix(read.csv("cts.csv", row.names = 1))

coldata <- read.csv("coldata.csv", row.names = 1)
coldata$batch <- as.factor(coldata$batch)
coldata <- coldata[, c("batch", "condition")]

colnames(cts) <- sub("X", "", colnames(cts))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# analysis with DESeq2 taking into account batch and condition
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ batch + condition)
dds

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Note on factor levels
dds$condition <- relevel(dds$condition, ref = "NT")

# Differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

# Extracting transformed values

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Effects of transformations on the variance

# this gives log2(n + 1)
ntd <- normTransform(dds)

library(vsn)
library(ggplot2)

par(mfrow = c(3,1))

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Heatmap of the count matrix

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("batch", "condition")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


# Heatmap of the sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$batch, vsd$condition, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$batch, vsd$condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col = colors)

# Principal component plot of the samples

plotPCA(vsd, intgroup=c("batch", "condition"))

pcaData <- plotPCA(vsd, intgroup=c("batch", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch, label = pcaData$condition)) +
  geom_point(size=3) +
  geom_text(color = "black", vjust = 0, nudge_y = 0.5) +
  ggtitle("PCA Vps37 project data") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  ggtitle("PCA Vps37 project data") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
