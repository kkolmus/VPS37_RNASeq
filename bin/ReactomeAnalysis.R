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
  library(ReactomePA)
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


#############################################
flog.debug("Analysis of signalling pathways")
#############################################

# load data

VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
genes.up = readRDS(file.path(proj.dir, data.dir, "genes.up.RDS"))
genes.down = readRDS(file.path(proj.dir, data.dir, "genes.down.RDS"))

# function

PathwayAnalyzer <- function(dataset, 
                            pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05,
                            minGSSize = 10, maxGSSize = 500, organism = "human") {
  enrichPathway(gene = dataset, 
                universe = VPS37_background[[2]], 
                organism = organism,
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff  = pvalueCutoff, 
                qvalueCutoff  = qvalueCutoff,
                minGSSize = minGSSize, 
                maxGSSize = maxGSSize, 
                readable = TRUE)
}

# for loop

Pathway.up <- setNames(vector("list", length(conditions)), conditions)

for (ds in names(genes.up)) {
  print(ds)
  Pathway.up[[ds]] <- PathwayAnalyzer(dataset = genes.up[[paste0(ds)]][[2]])
}

Pathway.up_final <- NULL

for (ds in names(Pathway.up)) {
  print(ds)
  if (nrow(as.data.frame(Pathway.up[[paste0(ds)]])) == 0) {
    next }
  else {
    temp = as.data.frame(Pathway.up[[paste0(ds)]])
    temp = mutate(temp, Type = factor(paste0(ds)))
    a = as.numeric(gsub("[0-9]+[:/:]", "", temp$GeneRatio))
    temp = mutate(temp, GeneRatio = Count/a)
    temp = filter(temp, Count >= 10) 
  }
  Pathway.up_final = rbind(Pathway.up_final, temp)
  Pathway.up_final = Pathway.up_final[!duplicated(Pathway.up_final), ]
}

# vector <- c(rep(0, 9), "siVPS37A")
# 
# Pathway.up_final$Type <- as.character(Pathway.up_final$Type)
# 
# class(Pathway.up_final$Type)
# 
# Pathway.up_final <- rbind(Pathway.up_final, vector)

list_of_pathways <- c("Signaling by Interleukins",
                      "Signaling by GPCR", 
                      "Signaling by Receptor Tyrosine Kinases",
                      "MAPK family signaling cascades",
                      "Interleukin-20 family signaling",
                      "Unfolded Protein Response (UPR)", 
                      #"PI5P, PP2A and IER3 Regulate PI3K/AKT Signaling",
                      "Interleukin-4 and 13 signaling",
                      "PI3K/AKT Signaling in Cancer")

Pathway.up_final <- filter(Pathway.up_final,
                           Pathway.up_final$Description %in% list_of_pathways)

rm(a, cutoff, ds, temp)

# plot pathways analysis

plotPathway <- ggplot(data = Pathway.up_final, 
                      mapping = aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits = c(10^(-10), 0.05), low = "red", high = "blue") +
  ylab(NULL) + xlab("Gene ratio") + 
  labs(colour = "p-value", size = "Gene count") + 
  # ggtitle(paste0("Signaling Pathways Enrichment\nUpregulated Signaling Pathways")) +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold", color = "black")) +
  theme(strip.text = element_text(colour = "black", size = 14)) +
  theme(axis.text.y = element_text(colour = "black", size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14))

p.pathway.up <- plotPathway + facet_grid(.~Type)
p.pathway.up


# save plot
tiff(filename = file.path(proj.dir, "figures", "p.pathways.up.tiff"), 
     width = 1000, height = 400, units = "px")
p.pathway.up
dev.off()


flog.debug("Session Info")
sessionInfo()
