rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(readxl)
  library(WriteXLS)
  
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  
  library(rlist)
  library(rlang)
  library(imputeTS)
  
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  library(biomaRt)
  library(ComplexHeatmap)
  fontsize <- 0.6
  library(VennDiagram)
  library(gridExtra)
  library(circlize)
  library(ReactomePA)
  library(ggrepel)
  library(seplyr)
  library(cowplot)
})


flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")
proj.dir = "~/Desktop/HT projects/VPS37_RNASeq"
input.dir = "input"
data.dir = "data"

DF <- read_excel(file.path(proj.dir, input.dir, "KK_Vps37_project_full_over100.xlsx"))



flog.debug("Set parameters for data analysis")
UP = 3/2         # minimal fold change to consider the gene to be upregulated
DOWN = 2/3       # minimal fold change to consider the gene to be downregulated
pval = 0.05      # p-value to consider change to be statistically significant
reads = 100      # minimal number of reads to consider the gene to be expressed

conditions <- c("siVPS37A#1", "siVPS37B#1", "siVPS37C#1", 
                "siVPS37AB#1", "siVPS37AC#1", "siVPS37BC#1", "siVPS37ABC#1")



flog.debug("Convert log2FC ===> FC")
DF <- mutate_at(DF, vars(-c(1,2, seq(from = 4, to = 34, by = 2), 35:44)), funs(2^(.)))



flog.debug("Rename columns")

new_colnames <- c()

for(c in colnames(DF)) {
  if(grepl("VPS37", c)) {
    temp = sub(".vs.", "#1.vs.", c) 
    print(temp)
    }
  else {
    temp = c
    print(temp)
  }
  new_colnames = c(new_colnames, temp)
}

rm(temo, c)

colnames(DF) <- new_colnames

rm(new_colnames)

# save data
saveRDS(DF, file.path(proj.dir, data.dir, "DF_step1.RDS"))


  
flog.debug("Drop genes with NA for adjusted p-value")

# read data
DF = readRDS(file.path(proj.dir, data.dir, "DF_step1.RDS"))

# function
drop_missing <- function(dataset, condition) {
  drop_na(data = dataset, 
          paste0(condition, ".vs.NT_adj.pval"), 
          paste0(condition, ".vs.siCtrl134_adj.pval"))
}

# empty list
dropped.missing <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (c in conditions) {
  print(c)
  dropped.missing[[c]] <- drop_missing(DF, c)
}

rm(c)

# save data
saveRDS(dropped.missing, file.path(proj.dir, data.dir, "DF_step2.RDS"))



flog.debug("Select genes with FC >= 1.5 and p-value < 0.05")

# read data
dropped.missing = readRDS(file.path(proj.dir, data.dir, "DF_step2.RDS"))

# function
filter.up <- function(dataset, condition) {
  condition_string_NT_FC = paste0(condition, ".vs.NT_FC")
  condition_string_NT_pval = paste0(condition, ".vs.NT_adj.pval")
  condition_string_CTRL_FC = paste0(condition, ".vs.siCtrl134_FC")
  condition_string_CTRL_pval = paste0(condition, ".vs.siCtrl134_adj.pval")
  condition_name_NT_FC <- rlang::sym(condition_string_NT_FC)
  condition_name_NT_pval <- rlang::sym(condition_string_NT_pval)
  condition_name_CTRL_FC <- rlang::sym(condition_string_CTRL_FC)
  condition_name_CTRL_pval <- rlang::sym(condition_string_CTRL_pval)
  filter(dataset, 
         UQ(condition_name_NT_FC) >= UP &
           UQ(condition_name_NT_pval) < pval &
           UQ(condition_name_CTRL_FC) >= UP &
           UQ(condition_name_CTRL_pval) < pval &
           readsSum >= reads)
}

# empty list
filtered.up <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(dropped.missing)) {
  print(ds)
  filtered.up[[ds]] <- filter.up(dataset = dropped.missing[[paste0(ds)]],
                                 condition = ds)
}

# save data
saveRDS(filtered.up, file.path(proj.dir, data.dir, "DF_step3_filtered.up.RDS"))



flog.debug("Select genes with FC <= 0.66 and p-value < 0.05")

# read data
dropped.missing = readRDS(file.path(proj.dir, data.dir, "DF_step2.RDS"))

# function
filter.down <- function(dataset, condition) {
  condition_string_NT_FC = paste0(condition, ".vs.NT_FC")
  condition_string_NT_pval = paste0(condition, ".vs.NT_adj.pval")
  condition_string_CTRL_FC = paste0(condition, ".vs.siCtrl134_FC")
  condition_string_CTRL_pval = paste0(condition, ".vs.siCtrl134_adj.pval")
  condition_name_NT_FC <- rlang::sym(condition_string_NT_FC)
  condition_name_NT_pval <- rlang::sym(condition_string_NT_pval)
  condition_name_CTRL_FC <- rlang::sym(condition_string_CTRL_FC)
  condition_name_CTRL_pval <- rlang::sym(condition_string_CTRL_pval)
  filter(dataset, 
         UQ(condition_name_NT_FC) <= DOWN &
           UQ(condition_name_NT_pval) < pval &
           UQ(condition_name_CTRL_FC) <= DOWN &
           UQ(condition_name_CTRL_pval) < pval &
           readsSum >= reads)
}

# empty list
filtered.down <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(dropped.missing)) {
  print(ds)
  filtered.down[[ds]] <- filter.down(dataset = dropped.missing[[paste0(ds)]],
                                     condition = ds)
}

# save data
saveRDS(filtered.down, file.path(proj.dir, data.dir, "DF_step3_filtered.down.RDS"))


flog.debug("Session Info")
sessionInfo()
