rm(list = ls())

suppressPackageStartupMessages({
  library(futile.logger)
  library(dplyr)
  library(tibble)
  library(purrr)
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
  library(ggrepel)
  library(seplyr)
})


flog.threshold(DEBUG)

flog.debug("Set working directory and load required data")
proj.dir = "~/Desktop/HT projects/VPS37_RNASeq"
input.dir = "input"
data.dir = "data"

conditions <- c("siVPS37A#1", "siVPS37B#1", "siVPS37C#1", 
                "siVPS37AB#1", "siVPS37AC#1", "siVPS37BC#1", "siVPS37ABC#1")

DF = readRDS(file.path(proj.dir, data.dir, "DF_step1.RDS"))
filtered.up = readRDS(file.path(proj.dir, data.dir, "DF_step3_filtered.up.RDS"))
filtered.down = readRDS(file.path(proj.dir, data.dir, "DF_step3_filtered.down.RDS"))
filtered = Map(rbind, filtered.up, filtered.down)

saveRDS(filtered, file.path(proj.dir, data.dir, "filtered.RDS"))


##############################
flog.debug("Convert gene IDs")
##############################

flog.debug("Selecting background genes for analysis")

VPS37_background <- bitr(DF$RefSeq, fromType = "REFSEQ", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# save data
saveRDS(VPS37_background, file.path(proj.dir, data.dir, "VPS37_background.RDS"))

flog.debug("Converting RefSeq to EntrezID for differentially expressed genes in each condition")

# function
convert.IDs <- function(dataset, OrgDb = "org.Hs.eg.db",
                        fromType = "REFSEQ", toType = "ENTREZID") {
  bitr(geneID = dataset, fromType = fromType, toType = toType,  OrgDb = OrgDb)
}

flog.debug("Differentially expressed genes")

# list
genes <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(filtered)) {
  print(ds)
  genes[[ds]] <- convert.IDs(dataset = filtered[[paste0(ds)]][[1]])
}

# save data
saveRDS(genes, file.path(proj.dir, data.dir, "genes.RDS"))





###########################################################################
flog.debug("GO analysis of biological processes for individual conditions")
###########################################################################

# read data
VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
genes = readRDS(file.path(proj.dir, data.dir, "genes.RDS"))

# function
GOanalyzer <- function(dataset, 
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05, qvalueCutoff  = 0.05,
                       minGSSize = 10, maxGSSize = 500) {
  enrichGO(gene = dataset, universe = VPS37_background[[2]], OrgDb = org.Hs.eg.db,
           ont = ont, pAdjustMethod = pAdjustMethod,
           pvalueCutoff  = pvalueCutoff, qvalueCutoff  = qvalueCutoff,
           minGSSize = minGSSize, maxGSSize = maxGSSize, readable = TRUE)
}

flog.debug("Biological processes for differentially expressed genes")

# list
BP <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(genes)) {
  print(ds)
  BP[[ds]] <- GOanalyzer(dataset = genes[[paste0(ds)]][[2]])
}

# save data
saveRDS(BP, file.path(proj.dir, data.dir, "BP.RDS"))


#########################################
flog.debug("Combine redundant processes")
#########################################

# read data

VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
BP = readRDS(file.path(proj.dir, data.dir, "BP.up.RDS"))

# function
GOsimplifier <- function(dataset, cutoff, by = "p.adjust") {
  clusterProfiler::simplify(x = dataset, cutoff = cutoff, by = by, select_fun = min)
}

flog.debug("Grouping of biological processes for differentially expressed genes")

# list
BP.sim.cutoff <- list()
BP.sim.condition <- NULL

# for loop
for (i in seq(from = 0.4, to = 0.70, by = 0.05)) {
  print(i)
  for (ds in names(BP)) {
    print(ds)
    if (nrow(as.data.frame(BP[[paste0(ds)]])) == 0) {
      next}
    else {
      temp = as.data.frame(GOsimplifier(dataset = BP[[paste0(ds)]], cutoff = i))
      temp = mutate(temp, Type = factor(paste0(ds)))
      a = as.numeric(gsub("[0-9]+[:/:]", "", temp$GeneRatio))
      temp = mutate(temp, GeneRatio = Count/a)
      temp = filter(temp, Count >= 10) }
    BP.sim.condition = rbind(BP.sim.condition, temp)
    BP.sim.condition = BP.sim.condition[!duplicated(BP.sim.condition), ]
  }
  BP.sim.cutoff = list.append(BP.sim.cutoff, BP.sim.condition)
}

names(BP.sim.cutoff) <- paste0("cutoff_", seq(from = 0.40, to = 0.70, by = 0.05))

# save data
saveRDS(BP.sim.cutoff, file.path(proj.dir, data.dir, "BP.simplified_list.RDS"))


############################
flog.debug("Plotting of BP")
############################

# load data

BP.sim.cutoff = readRDS(file.path(proj.dir, data.dir, "BP.simplified_list.RDS"))

# function

plotting <- function(data, alteration) {
  
  tempABC = filter(data, Type == "siVPS37ABC#1")
  tempABC = top_n(x = tempABC, n = 15, wt = Count)
  tempABC = tempABC[, 1]
  
  input = filter(data, data[["ID"]] %in% tempABC)
  
  plotBP <- ggplot(data = input,
                   mapping = aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) + 
    scale_colour_gradient(limits = c(10^(-10), 0.05), low = "red", high = "blue") +
    ylab(NULL) + xlab("Gene ratio") + 
    labs(colour = "p-value", size = "Gene count") + 
    # ggtitle(paste0("Gene Ontology Enrichment\n", alteration, "regulated Biological Processes")) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold", color = "black")) +
    theme(strip.text = element_text(colour = "black", size = 14)) +
    theme(axis.text.y = element_text(colour = "black", size = 14)) +
    theme(axis.text.x = element_text(colour = "black", size = 14))
  
  plotBP + facet_grid(.~Type)
}

# print plots for different cutoff values

for (cutoff in names(BP.sim.cutoff)) {
  print(cutoff)
  print(plotting(data = BP.sim.cutoff[[paste0(cutoff)]]))
}

# print plots for the best cutoff value

p  <- plotting(data = BP.sim.cutoff[["cutoff_0.6"]])
p

# save plot
tiff(filename = file.path(proj.dir, "figures", "p.up.tiff"), 
     width = 1000, height = 400, units = "px")
p
dev.off()


flog.debug("Session Info")
sessionInfo()
