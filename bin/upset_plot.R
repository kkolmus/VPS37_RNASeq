# install.packages("UpSetR")

rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(tidyverse)
  library(UpSetR)
})



flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")

proj.dir = "~/Desktop/HT projects/VPS37_RNASeq"
input.dir = "input"
data.dir = "data"

conditions <- c("siVPS37A#1", "siVPS37B#1", "siVPS37C#1", 
                "siVPS37AB#1", "siVPS37AC#1", "siVPS37BC#1", "siVPS37ABC#1")

filter_final <- readRDS(file.path(proj.dir, data.dir, "filtered.final.RDS"))

gene.list <- setNames(vector("list", length(conditions)), conditions)

for (i in 1:length(filter_final)) {
  temp_df = filter_final[[i]][["Overlap"]]
  temp_name = names(filter_final[i])
  gene.list[[temp_name]] = temp_df
}

upset(fromList(gene.list),  
      nsets = 7, 
      nintersects = NA, 
      order.by = "freq",
      #keep.order = TRUE,
      mainbar.y.label = "\n \n \nNumber of genes", # "Intersection size \nof differentially expressed genes", 
      sets.x.label = "", # "Number of differentially \nexpressed genes"
      text.scale = c(7, 7, 7, 0, 7, 4),
      mainbar.y.max = 600,
      point.size = 5,
      #empty.intersections = "on",
      set_size.show = TRUE,
      number.angles = 0
      )
