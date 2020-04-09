rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(rlist)
  library(rlang)
  library(imputeTS)
  library(VennDiagram)
  library(gridExtra)
})


flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")

proj.dir = "~/Desktop/HT projects/VPS37_RNASeq"
input.dir = "input"
data.dir = "data"

conditions <- c("siVPS37A#1", "siVPS37B#1", "siVPS37C#1", 
                "siVPS37AB#1", "siVPS37AC#1", "siVPS37BC#1", "siVPS37ABC#1")

UP = 3/2         # minimal fold change to consider the gene to be upregulated
DOWN = 2/3       # minimal fold change to consider the gene to be downregulated
pval = 0.05      # p-value to consider change to be statistically significant
reads = 100      # minimal number of reads to consider the gene to be expressed

dropped.missing = readRDS(file.path(proj.dir, data.dir, "DF_step2.RDS"))


########################### 
flog.debug("Venn diagrams")
###########################

# Upregulated genes

# functions

filter.NT.function <- function(dataset, condition) {
  condition_string_NT_FC = paste0(condition, ".vs.NT_FC")
  condition_string_NT_pval = paste0(condition, ".vs.NT_adj.pval")
  condition_name_NT_FC <- rlang::sym(condition_string_NT_FC)
  condition_name_NT_pval <- rlang::sym(condition_string_NT_pval)
  filter(dataset, 
         UQ(condition_name_NT_FC) >= UP |
           UQ(condition_name_NT_FC) <= DOWN &
           UQ(condition_name_NT_pval) < pval &
           readsSum >= reads)
}

filter.CTRL.function <- function(dataset, condition) {
  condition_string_CTRL_FC = paste0(condition, ".vs.siCtrl134_FC")
  condition_string_CTRL_pval = paste0(condition, ".vs.siCtrl134_adj.pval")
  condition_name_CTRL_FC <- rlang::sym(condition_string_CTRL_FC)
  condition_name_CTRL_pval <- rlang::sym(condition_string_CTRL_pval)
  filter(dataset, 
         UQ(condition_name_CTRL_FC) >= UP |
           UQ(condition_name_CTRL_FC) <= DOWN &
           UQ(condition_name_CTRL_pval) < pval &
           readsSum >= reads)
}


# for loog

filtered.final <- NULL

for (ds in names(dropped.missing)) {
  print(ds)
  filtered.final[[ds]][["siCTRL"]] <- filter.CTRL.function(
    dataset = dropped.missing[[paste0(ds)]], condition = ds)[,1]
  filtered.final[[ds]][["NT"]] <- filter.NT.function(
    dataset = dropped.missing[[paste0(ds)]], condition = ds)[,1]
  filtered.final[[ds]][["Overlap"]] <- intersect(
    filtered.final[[ds]][["siCTRL"]][["RefSeq"]],
    filtered.final[[ds]][["NT"]][["RefSeq"]])
}

draw.2.Venn <- function (area1, area2, intersection, condition) {
  draw.pairwise.venn(area1 = length(area1),
                     area2 = length(area2),
                     cross.area = length(intersection),
                     # c(paste0(condition, "/NT"), paste0(condition, "/siCTRL")),
                     fill = c("dodgerblue1", "deeppink2"),
                     fontfamily = "Arial", cat.fontfamily = "Arial",
                     cex = 0, #
                     # cat.cex = 4, #
                     # cat.pos = c(30, 215), #
                     # cat.dist = .05, #
                     scaled = TRUE)
}

for (i in names(filtered.final)) {
  print(i)
  grid.newpage()
  print(draw.2.Venn(area1 = filtered.final[[paste0(i)]][["siCTRL"]][["RefSeq"]],
                    area2 = filtered.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.final[[paste0(i)]][["Overlap"]],
                    condition = i))
  print(draw.2.Venn(area1 = filtered.final[[paste0(i)]][["siCTRL"]][["RefSeq"]],
                    area2 = filtered.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.final[[paste0(i)]][["Overlap"]],
                    condition = i))
  print(draw.2.Venn(area1 = filtered.final[[paste0(i)]][["siCTRL"]][["RefSeq"]],
                    area2 = filtered.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.final[[paste0(i)]][["Overlap"]],
                    condition = i))
  print(draw.2.Venn(area1 = filtered.final[[paste0(i)]][["siCTRL"]][["RefSeq"]],
                    area2 = filtered.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.final[[paste0(i)]][["Overlap"]],
                    condition = i))
  print(draw.2.Venn(area1 = filtered.final[[paste0(i)]][["siCTRL"]][["RefSeq"]],
                    area2 = filtered.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.final[[paste0(i)]][["Overlap"]],
                    condition = i))
}

# save data
saveRDS(filtered.final, file.path(proj.dir, data.dir, "filtered.final.RDS"))


flog.debug("Session Info")
sessionInfo()