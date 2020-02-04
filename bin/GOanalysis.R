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

flog.debug("Upregulated genes")

# list
genes.up <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(filtered.up)) {
  print(ds)
  genes.up[[ds]] <- convert.IDs(dataset = filtered.up[[paste0(ds)]][[1]])
}

# save data
saveRDS(genes.up, file.path(proj.dir, data.dir, "genes.up.RDS"))

flog.debug("Downregulated genes")

# list
genes.down <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(filtered.up)) {
  print(ds)
  genes.down[[ds]] <- convert.IDs(dataset = filtered.down[[paste0(ds)]][[1]])
}

# save data
saveRDS(genes.down, file.path(proj.dir, data.dir, "genes.down.RDS"))



###########################################################################
flog.debug("GO analysis of biological processes for individual conditions")
###########################################################################

# read data
VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
genes.up = readRDS(file.path(proj.dir, data.dir, "genes.up.RDS"))
genes.down = readRDS(file.path(proj.dir, data.dir, "genes.down.RDS"))

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

flog.debug("Biological processes for transcriptionally upregulated genes")

# list
BP.up <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(genes.up)) {
  print(ds)
  BP.up[[ds]] <- GOanalyzer(dataset = genes.up[[paste0(ds)]][[2]])
}

# save data
saveRDS(BP.up, file.path(proj.dir, data.dir, "BP.up.RDS"))

flog.debug("Biological processes for transcriptionally downregulated genes")

# list
BP.down <- setNames(vector("list", length(conditions)), conditions)

# for loop
for (ds in names(genes.up)) {
  print(ds)
  BP.down[[ds]] <- GOanalyzer(dataset = genes.down[[paste0(ds)]][[2]])
}

# save data
saveRDS(BP.down, file.path(proj.dir, data.dir, "BP.down.RDS"))



#########################################
flog.debug("Combine redundant processes")
#########################################

# read data

VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
BP.up = readRDS(file.path(proj.dir, data.dir, "BP.up.RDS"))
BP.down = readRDS(file.path(proj.dir, data.dir, "BP.down.RDS"))

# function
GOsimplifier <- function(dataset, cutoff, by = "p.adjust") {
  clusterProfiler::simplify(x = dataset, cutoff = cutoff, by = by, select_fun = min)
}

flog.debug("Grouping of biological processes for transcriptionally upregulated genes")

# list
BP.up.sim.cutoff <- list()
BP.up.sim.condition <- NULL

# for loop
for (i in seq(from = 0.4, to = 0.70, by = 0.05)) {
  print(i)
  for (ds in names(BP.up)) {
    print(ds)
    if (nrow(as.data.frame(BP.up[[paste0(ds)]])) == 0) {
      next}
    else {
      temp = as.data.frame(GOsimplifier(dataset = BP.up[[paste0(ds)]], cutoff = i))
      temp = mutate(temp, Type = factor(paste0(ds)))
      a = as.numeric(gsub("[0-9]+[:/:]", "", temp$GeneRatio))
      temp = mutate(temp, GeneRatio = Count/a)
      temp = filter(temp, Count >= 10) }
    BP.up.sim.condition = rbind(BP.up.sim.condition, temp)
    BP.up.sim.condition = BP.up.sim.condition[!duplicated(BP.up.sim.condition), ]
  }
  BP.up.sim.cutoff = list.append(BP.up.sim.cutoff, BP.up.sim.condition)
}

names(BP.up.sim.cutoff) <- paste0("cutoff_", seq(from = 0.40, to = 0.70, by = 0.05))

# save data
saveRDS(BP.up.sim.cutoff, file.path(proj.dir, data.dir, "BP.up.simplified_list.RDS"))


flog.debug("Grouping of biological processes for transcriptionally downregulated genes")

# list
BP.down.sim.cutoff <- list()
BP.down.sim.condition <- NULL

# for loop
for (i in seq(from = 0.4, to = 0.70, by = 0.05)) {
  print(i)
  for (ds in names(BP.down)) {
    print(ds)
    if (nrow(as.data.frame(BP.down[[paste0(ds)]])) == 0) {
      next}
    else {
      temp = as.data.frame(GOsimplifier(dataset = BP.down[[paste0(ds)]], cutoff = i))
      temp = mutate(temp, Type = factor(paste0(ds)))
      a = as.numeric(gsub("[0-9]+[:/:]", "", temp$GeneRatio))
      temp = mutate(temp, GeneRatio = Count/a)
      temp = filter(temp, Count >= 10) }
    BP.down.sim.condition = rbind(BP.down.sim.condition, temp)
    BP.down.sim.condition = BP.down.sim.condition[!duplicated(BP.down.sim.condition), ]
  }
  BP.down.sim.cutoff = list.append(BP.down.sim.cutoff, BP.down.sim.condition)
}

names(BP.down.sim.cutoff) <- paste0("cutoff_", seq(from = 0.40, to = 0.70, by = 0.05))

# save data
saveRDS(BP.down.sim.cutoff, file.path(proj.dir, data.dir, "BP.down.simplified_list.RDS"))

rm(temp, a, i, ds, BP.down.sim.condition, BP.up.sim.condition)



############################
flog.debug("Plotting of BP")
############################

# load data

BP.up.sim.cutoff = readRDS(file.path(proj.dir, data.dir, "BP.up.simplified_list.RDS"))
BP.down.sim.cutoff = readRDS(file.path(proj.dir, data.dir, "BP.down.simplified_list.RDS"))

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

for (cutoff in names(BP.up.sim.cutoff)) {
  print(cutoff)
  print(plotting(data = BP.up.sim.cutoff[[paste0(cutoff)]], alteration = "Up"))
}

# print plots for the best cutoff value

p.up  <- plotting(data = BP.up.sim.cutoff[["cutoff_0.6"]], alteration = "Up")
p.up
p.down <- plotting(data = BP.down.sim.cutoff[["cutoff_0.6"]], alteration = "Down")
p.down

# save plot
tiff(filename = file.path(proj.dir, "figures", "p.up.tiff"), 
     width = 1000, height = 400, units = "px")
p.up
dev.off()


flog.debug("Session Info")
sessionInfo()
