# cleaning
rm(list = ls())

# Normalisation vs. NT
setwd("~/Desktop/HT projects/VPS37_RNASeq/NT")

files <- list.files(path="~/Desktop/HT projects/VPS37_RNASeq/NT")
genes <- read.table(files[1], row.names = 1, header=TRUE, sep="\t")[,0] # gene names
df <- do.call(cbind,
              lapply(files, function(fn) 
                read.table(fn, header=TRUE, sep="\t")[,c(2,6)]))
prefix <- gsub(".txt", "", files)
suffix <- c("FC", "adj.pval")
DF_NT <- cbind(genes,df)

setwd("~/Desktop/HT projects/VPS37_RNASeq")
write.table(DF_NT, "DF_NT.txt", quote = F, sep = "\t", row.names = T)

# cleaning
rm(list = ls())

# Normalisation vs. siCtrl
setwd("~/Desktop/HT projects/VPS37_RNASeq/siCTRL")

files <- list.files(path="~/Desktop/HT projects/VPS37_RNASeq/siCTRL")
genes <- read.table(files[1], row.names = 1, header=TRUE, sep="\t")[,0] # gene names
df <- do.call(cbind,
              lapply(files, function(fn) 
                read.table(fn, header=TRUE, sep="\t")[,c(2,6)]))
prefix <- gsub(".txt", "", files)
suffix <- c("FC", "adj.pval")
DF_siCTRL <- cbind(genes,df)

setwd("~/Desktop/HT projects/VPS37_RNASeq")
write.table(DF_siCTRL, "DF_siCTRL.txt", quote = F, sep = "\t", row.names = T)

# cleaning
rm(list = ls())

# edit datasheets

setwd("~/Desktop/HT projects/VPS37_RNASeq")

library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(imputeTS)

# merging data

DF_NT <- read_excel("DF_NT.xlsx")
DF_siCTRL <- read_excel("DF_siCTRL.xlsx")

DF <- merge(DF_NT, DF_siCTRL, by = "RefSeq", all = TRUE)

# add gene symbols to RefSeq

DF_Symbol <- bitr(DF$RefSeq, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

new_col_names2 <- c("RefSeq", "Symbol")

colnames(DF_Symbol) <- new_col_names2

DF <- merge(DF, DF_Symbol, by = "RefSeq", all = TRUE)

# change class of columns

DF <- mutate_at(DF, vars(-c(RefSeq, Symbol)), funs(as.numeric(.)))

# change log2FC to FC

DF2 <- mutate_at(DF, vars(-c(RefSeq, Symbol,
                             siCTRL.vs.NT_adj.pval,
                             siVPS37A_1.vs.NT_adj.pval, 
                             siVPS37B_1.vs.NT_adj.pval, 
                             siVPS37C_1.vs.NT_adj.pval, 
                             siVPS37AB_1.vs.NT_adj.pval,
                             siVPS37AC_1.vs.NT_adj.pval,
                             siVPS37BC_1.vs.NT_adj.pval,
                             siVPS37ABC_1.vs.NT_adj.pval,
                             NT.vs.siCTRL_adj.pval,
                             siVPS37A_1.vs.siCTRL_adj.pval,
                             siVPS37B_1.vs.siCTRL_adj.pval,
                             siVPS37C_1.vs.siCTRL_adj.pval,
                             siVPS37AB_1.vs.siCTRL_adj.pval,
                             siVPS37AC_1.vs.siCTRL_adj.pval,
                             siVPS37BC_1.vs.siCTRL_adj.pval,
                             siVPS37ABC_1.vs.siCTRL_adj.pval)), 
                 funs(2^(.)))

# round up numbers

DF2 <- mutate_at(DF2, vars(-c(RefSeq, Symbol,
                              siCTRL.vs.NT_adj.pval,
                              siVPS37A_1.vs.NT_adj.pval, 
                              siVPS37B_1.vs.NT_adj.pval, 
                              siVPS37C_1.vs.NT_adj.pval, 
                              siVPS37AB_1.vs.NT_adj.pval,
                              siVPS37AC_1.vs.NT_adj.pval,
                              siVPS37BC_1.vs.NT_adj.pval,
                              siVPS37ABC_1.vs.NT_adj.pval,
                              NT.vs.siCTRL_adj.pval,
                              siVPS37A_1.vs.siCTRL_adj.pval,
                              siVPS37B_1.vs.siCTRL_adj.pval,
                              siVPS37C_1.vs.siCTRL_adj.pval,
                              siVPS37AB_1.vs.siCTRL_adj.pval,
                              siVPS37AC_1.vs.siCTRL_adj.pval,
                              siVPS37BC_1.vs.siCTRL_adj.pval,
                              siVPS37ABC_1.vs.siCTRL_adj.pval)), 
                 funs(round(., 3)))

# add avaeraged number of reads per gene

reads <- read.csv("cts.csv", header = TRUE)
reads <- mutate(reads, 
                NT_reads = (X1_NT + X2_NT + X3_NT)/3,
                siCTRL_reads = (X1_siCTRL + X2_siCTRL + X3_siCTRL)/3,
                siVPS37ABC_1_reads = (X1_siVPS37ABC_1 + X2_siVPS37ABC_1 + X3_siVPS37ABC_1)/3,
                siVPS37AB_1_reads = (X1_siVPS37AB_1 + X2_siVPS37AB_1 + X3_siVPS37AB_1)/3,
                siVPS37AC_1_reads = (X1_siVPS37AC_1 + X2_siVPS37AC_1 + X3_siVPS37AC_1)/3,
                siVPS37BC_1_reads = (X1_siVPS37BC_1 + X2_siVPS37BC_1 + X3_siVPS37BC_1)/3,
                siVPS37A_1_reads = (X1_siVPS37A_1 + X2_siVPS37A_1 + X3_siVPS37A_1)/3,
                siVPS37B_1_reads = (X1_siVPS37B_1 + X2_siVPS37B_1 + X3_siVPS37B_1)/3,
                siVPS37C_1_reads = (X1_siVPS37C_1 + X2_siVPS37C_1 + X3_siVPS37C_1)/3)
reads <- mutate_at(reads, vars(-c(RefSeq)), funs(round(., digits = 3)))

DF2 <- left_join(DF2, reads, by = "RefSeq")

# further filtering

library(readxl)
library(dplyr)
library(WriteXLS)

DF2 <- mutate(DF2, readsSum = rowSums(DF2[, c(35:44)]))

DF_filtered <- filter(DF2, readsSum >= 10)

WriteXLS(DF_filtered, "KK_Vps37_project_full.xlsx")

DF_filtered <- filter(DF2, readsSum >= 100)

WriteXLS(DF_filtered, "KK_Vps37_project_full_over100.xlsx")

