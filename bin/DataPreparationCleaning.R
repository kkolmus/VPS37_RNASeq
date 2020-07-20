# Data frame preparation for Tsg101/Vps28 project

# cleaning
rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(DESeq2)
})

# setting up work space for analysis
setwd("~/Desktop/HT projects/VPS37_RNASeq")

# prepare count matrix for analysis
DF <- read_excel("rawReads.Vps37.xlsx")

# remove statistics for columns and hgnc column
DF <- DF[-c(1:10,20813:20817),-c(2)]

# fill up missing values
DF[is.na(DF)] <- 0

# remove duplicates
DF <- DF[! duplicated(DF$RefSeq),] # removes 0 indices

# exporting data frame
write.csv(DF, "cts.csv", row.names = FALSE)

# loading data on experimental design
coldata <- read_excel("coldata.xlsx")

# exporting data frame
write.csv(coldata, "coldata.csv", row.names = FALSE)

# loading data frames for DESeq2 analysis
cts <- as.matrix(read.csv("cts.csv", row.names = 1))

coldata <- read.csv("coldata.csv", row.names = 1)
coldata$batch <- as.factor(coldata$batch)
coldata <- coldata[, c("batch", "condition")]

# confirming coherence

colnames(cts) <- sub("X", "", colnames(cts))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))