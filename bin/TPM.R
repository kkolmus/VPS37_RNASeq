### This script aims at preparing a matrix of normalized counts for the submission of the VPS37 data to GEO.

### LIBS ###

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(WriteXLS)
  library(biomaRt)
})

### FUNCTION ###

listMarts()
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

TPM_calculator <- function(dataset = DF) {
  IDs <- dataset$RefSeq
  length = getBM(attributes = c("refseq_mrna", "transcript_length"), 
                 filters = "refseq_mrna", values = IDs, mart = ensembl)
  length = length[! duplicated(length$refseq_mrna),]
  length = dplyr::select(length, RefSeq = refseq_mrna, transcript_length)
  df_temp = merge(dataset, length, by = c("RefSeq"))
  df_temp = df_temp[, c(1, 2, ncol(df_temp), 3:(ncol(df_temp)-1))]
  # normalize for gene length, TPM1 = reads per kilobase = RPK
  TPM1 <- mutate_at(df_temp, vars(-c(RefSeq, Symbol)), funs(./df_temp$transcript_length))[, -3]
  # normalize for sequencing depth, TPM2 = scaling factors
  TPM2 <- unname(unlist(apply(TPM1[, c(3:ncol(TPM1))], 2, sum)/10^(6)))
  # transcripts per million
  TPM_temp <- sweep(TPM1[,-c(1,2)], MARGIN = 2, TPM2, FUN = "/")
  TPM_names <- TPM1[, c(1,2)]
  TPM <<- cbind(TPM_names, TPM_temp) # TPM
  # scaled samples
  Av <- apply(TPM[, c(3:ncol(TPM))], 1, mean)
  StDev <- apply(TPM[, c(3:ncol(TPM))], 1, sd)
  TPM_HM <<- mutate_at(TPM, vars(-c(RefSeq, Symbol)), funs(((.)-Av)/StDev)) # scaled
}


### MAIN ###

dir <- "~/Desktop/HT projects/RNASeq Vps37/Complete analysis - only strong"

DF.TPM_normalization <- readRDS("~/Desktop/HT projects/VPS37_RNASeq/data/DF.TPM_normalization.RDS")
#DF.TPM_normalization <- DF.TPM_normalization[, c(1,2)]

DF.raw_reads <- readxl::read_xlsx(file.path(dir, "rawReads.Vps37.xlsx"))
DF.raw_reads <- dplyr::rename(DF.raw_reads, Symbol = symbol)

DF <- inner_join(DF.raw_reads, DF.TPM_normalization, by = c("RefSeq", "Symbol"))
DF <- DF[, -c(10,11, 21,22, 32,33)]
DF <- DF[, -c(30:38)]

old.colnames <- colnames(DF)
new.colnames <- c("RefSeq", "Symbol",
                  "1_siVPS37ABC#1", "1_siVPS37AB#1", "1_siVPS37AC#1", "1_siVPS37BC#1",
                  "1_siVPS37A#1", "1_siVPS37B#1", "1_siVPS37C#1", "1_siCTRL#1", "1_NT",
                  "2_siVPS37ABC#1", "2_siVPS37AB#1", "2_siVPS37AC#1", "2_siVPS37BC#1",
                  "2_siVPS37A#1", "2_siVPS37B#1", "2_siVPS37C#1", "2_siCTRL#1", "2_NT",
                  "3_siVPS37ABC#1", "3_siVPS37AB#1", "3_siVPS37AC#1", "3_siVPS37BC#1",
                  "3_siVPS37A#1", "3_siVPS37B#1", "3_siVPS37C#1", "3_siCTRL#1", "3_NT" )

colnames(DF) <- new.colnames

### Prepare matrix


TPM_calculator(DF)

saveRDS(TPM, "~/Desktop/HT projects/VPS37_RNASeq/data/DFsamples.TPM_normalization.RDS")