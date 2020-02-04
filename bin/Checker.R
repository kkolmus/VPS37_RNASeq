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
})


flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")
proj.dir = "~/Desktop/HT projects/VPS37_RNASeq"
input.dir = "input"
data.dir = "data"

conditions <- c("siVPS37A#1", "siVPS37B#1", "siVPS37C#1", 
                "siVPS37AB#1", "siVPS37AC#1", "siVPS37BC#1", "siVPS37ABC#1")

DF = readRDS(file.path(proj.dir, data.dir, "DF_step1.RDS"))
filtered.up.final = readRDS(file.path(proj.dir, data.dir, "filtered.up.final.RDS"))
filtered.down.final = readRDS(file.path(proj.dir, data.dir, "filtered.down.final.RDS"))


##############################################################################################
flog.debug("Prepare a list of dataframe with column pointing out significance of alterations")
##############################################################################################

checker <- function(dataset, condition, normalization, UP = 1.5, DOWN = 0.5, pval = 0.05) {
  condition_FC = paste0("si", condition, ".vs.", normalization, "_FC")
  condition_pval = paste0("si", condition, ".vs.", normalization, "_adj.pval")
  dataset = as.data.frame(mutate(dataset, threshold = with(dataset, 
                                                           ifelse(dataset[, condition_FC] >= UP &
                                                                    dataset[, condition_pval] < pval, "Upregulated",
                                                                  ifelse(dataset[, condition_FC] <= DOWN &
                                                                           dataset[, condition_pval] < pval,
                                                                         "Downregulated", "Not significant")))))
}

DF_checked <- setNames(vector("list", length(conditions)), conditions)

for (c in conditions) {
  DF_checked[[c]][["NT"]] <- checker(dataset = DF, condition = c, normalization = "NT")
  DF_checked[[c]][["siCtrl134"]] <- checker(dataset = DF, condition = c, normalization = "siCtrl134")
}

##################################################################################
flog.debug("Provide a list of genes you wish to visualize using the Volcano Plot")
##################################################################################

genes.input <- filter(DF, Symbol %in% c("insert genes of interest"))

for (i in names(filtered.up.final)) {
  print(i)
  grid.newpage()
  print(draw.2.Venn(area1 = filtered.up.final[[paste0(i)]][["CTRL"]][["RefSeq"]],
                    area2 = filtered.up.final[[paste0(i)]][["NT"]][["RefSeq"]],
                    intersection = filtered.up.final[[paste0(i)]][["Overlap"]],
                    condition = i))
}

# function to generate volcano plots

volcano.plot.generator <- function (dataset = DF_checked, condition, input = genes.input) {
  genes.input <- filter(data = DF_checked[[condition]][["NT"]], 
                        Symbol %in% c("genes of interest"))
  a = paste0("si", condition, ".vs.NT_FC")
  b = paste0("si", condition, ".vs.NT_adj.pval")
  vp1 <- ggplot(data = DF_checked[[condition]][["NT"]], 
                mapping = aes(x = log2(DF_checked[[condition]][["NT"]][[a]]),
                              y = -log10(DF_checked[[condition]][["NT"]][[b]]), 
                              colour = threshold)) +
    scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
    geom_point(alpha = 0.4, size = 1.0) + xlim(c(-4, 4)) + ylim(c(0, 50)) + 
    labs(color = "Expression pattern") +
    geom_point(data = genes.input, colour = "black") +
    ggtitle(paste0("si", condition, "/NT")) + 
    xlab("log2FoldChange") + ylab("-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(vp1 + 
          geom_text_repel(data = genes.input, mapping = aes(label = genes.input[["Symbol"]]), 
                          colour = "black", size = 4))
}

# for loop 

for (i in names(DF_checked)) {
  print(i)
  print(volcano.plot.generator(condition = i))
}


flog.debug("Session Info")
sessionInfo()