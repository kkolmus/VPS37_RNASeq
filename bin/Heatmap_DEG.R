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
  library(circlize)
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

DF = readRDS(file.path(proj.dir, data.dir, "DF_step1.RDS"))
VPS37_background = readRDS(file.path(proj.dir, data.dir, "VPS37_background.RDS"))
genes = readRDS(file.path(proj.dir, data.dir, "genes.RDS"))
BP.sim.cutoff = readRDS(file.path(proj.dir, data.dir, "BP.simplified_list.RDS"))


###################################################################
flog.debug("Dataframe with differentially expressed genes in the given condition")
###################################################################

sfx <- conditions
res <- genes[[1]]

for(i in head(seq_along(genes), -1)) {
  res <- merge(res, genes[[i+1]], all = TRUE, 
                  suffixes = sfx[i:(i+1)], by = "REFSEQ")
}

res <- res[,1]
res <- as.vector(res)


#########################################
flog.debug("Transcripts Per Million TPM")
#########################################

samples <- DF[, c(1,2, grep("_reads", colnames(DF)))]
names(samples) <- gsub(x = names(samples), pattern = "_reads", replacement = "")
colnames(samples) <- paste(colnames(samples), "#1", sep = "")
samples <- dplyr::rename(samples, RefSeq = `RefSeq#1`)
samples <- dplyr::rename(samples, Symbol = `Symbol#1`)
samples <- dplyr::rename(samples, NT = `NT#1`)
samples <- dplyr::rename(samples, `siCTRL#1` = `siCtrl134#1`)
samples <- samples[!grepl("NR_", samples$RefSeq),]

#sum(is.na(samples$RefSeq))                 # check if there're missing values
#samples[duplicated(samples$RefSeq),]       # check if there're duplicated rows

listMarts()
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# listDatasets(ensembl)

flog.debug("Prepate TPM dataframe")

TPM_calculator <- function(dataset = samples) {
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
  TPM <<- cbind(TPM_names, TPM_temp)
  # scaled samples
  Av <- apply(TPM[, c(3:ncol(TPM))], 1, mean)
  StDev <- apply(TPM[, c(3:ncol(TPM))], 1, sd)
  TPM_HM <<- mutate_at(TPM, vars(-c(RefSeq, Symbol)), funs(((.)-Av)/StDev))
}

TPM_calculator(samples)

# save data

saveRDS(TPM_HM, file.path(proj.dir, data.dir, "DF.TPM_normalization_UPandDOWNgenes.RDS"))


flog.debug("Heatmap for differentially regulated genes")

# laod data
TPM_HM = readRDS(file.path(proj.dir, data.dir, "DF.TPM_normalization_UPandDOWNgenes.RDS"))

# filter data to select upregulated genes
HeatmapData <- filter(TPM_HM, RefSeq %in% res)
HeatmapData_matrix <- as.matrix(HeatmapData[ , c(3:11)])
HeatmapData_matrix <- t(HeatmapData_matrix)

# draw heatmap
Heatmap(HeatmapData_matrix,
        cluster_columns = TRUE,
        column_names_side = "top", column_dend_side = "top",
        row_names_side = "left", row_dend_side = "left",
        row_names_gp = gpar(fontsize = 14), # fontface = "bold"
        row_dend_width = unit(2, "cm"),
        clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
        name = "z-score",
        column_title_gp = gpar(fontsize = 14, fontfamily = "Helvenica"), # fontface = "bold",
        row_title_rot = 0,
        col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 40)) # fontface = "bold")


flog.debug("Heatmap for downregulated genes")

HeatmapData <- filter(TPM_HM, RefSeq %in% res.down)
HeatmapData_matrix <- as.matrix(HeatmapData[ , c(3:11)])
HeatmapData_matrix <- t(HeatmapData_matrix)

Heatmap(HeatmapData_matrix,
        cluster_columns = TRUE,
        column_names_side = "top", column_dend_side = "top",
        row_names_side = "left", row_dend_side = "left",
        row_names_gp = gpar(fontsize = 14, fontface = "bold"),
        row_dend_width = unit(2, "cm"),
        clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
        name = "z-score",
        column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"), 
        row_title_rot = 0,
        col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
        row_title_gp = gpar(fontsize = 40, fontface = "bold")
        )


###############################
flog.debug("Heatmap generator")
###############################

# function

heatmap.generator <- function(
  dataset = DF, # dataframe with a list of genes, FC and adj. p-values 
  TPM_HM, # dataframe with a list of genes and read counts presented as TPM
  list.of.BP = BP.sim.cutoff, # dataframe with biological process
  GOid = "GO:0006954", # gene ontology identified e.g. "GO:0006954"
  condition = "siVPS37ABC", # condition for which biological processes will be presented
  cutoff = 0.55 # numeric input to access the list of simplified dataframe with biological processes
) {
  process = filter(list.of.BP[[paste0("cutoff_", cutoff)]], 
                   ID %in% GOid & Type == paste0(condition))
  process.list = process$geneID
  process.string = unlist(strsplit(process.list, split="/"))
  
  process.heatmap <- filter(dataset, Symbol %in% process.string)
  process.heatmap <- filter(TPM_HM, RefSeq %in% process.heatmap[[1]])
  process.input <- process.heatmap[, c(3:ncol(process.heatmap))]
  rownames(process.input) <- process.heatmap[,2]
  process.input <- t(process.input)
  
  Heatmap(process.input,
          cluster_columns = TRUE,
          column_names_side = "top", column_dend_side = "top",
          row_names_side = "left", row_dend_side = "left",
          row_names_gp = gpar(fontsize = 14),
          row_dend_width = unit(2, "cm"),
          clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
          name = "z-score",
          column_title = paste0(process$Description, " after si", condition),
          #row_title = "Silencing conditions",
          column_title_gp = gpar(fontsize = 14, fontfamily = "Helvenica"), # fontface = "bold",
          row_title_rot = 0,
          col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
          rect_gp = gpar(col = "black"),
          row_title_gp = gpar(fontsize = 40)) # fontface = "bold"
}


flog.debug("inflammatory response")

heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0006954", condition = "siVPS37ABC#1", cutoff = 0.6)

heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0006954", condition = "siVPS37AB#1", cutoff = 0.6)


flog.debug("regulation of growth")

heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0040008", condition = "siVPS37AB#1", cutoff = 0.6)

heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0040008", condition = "siVPS37AB#1", cutoff = 0.6)



# modification of function

total.process.string <- c()

heatmap.generator <- function(
  dataset = DF, # dataframe with a list of genes, FC and adj. p-values 
  TPM_HM, # dataframe with a list of genes and read counts presented as TPM
  list.of.BP = BP.sim.cutoff, # dataframe with biological process
  GOid = "GO:0006954", # gene ontology identified e.g. "GO:0006954"
  condition = c("siVPS37A#1", "siVPS37AB#1", "siVPS37AC#1", "siVPS37ABC#1"), # condition for which biological processes will be presented
  cutoff = 0.6 # numeric input to access the list of simplified dataframe with biological processes
) 
{
  for(c in condition) {
    print(c)
    process <- filter(list.of.BP[[paste0("cutoff_", cutoff)]], 
                      ID %in% GOid, Type == paste0(c))
    process.list = process$geneID
    process.string = unlist(strsplit(process.list, split="/"))
    total.process.string <<- c(total.process.string, process.string)
  }
  
  process.heatmap <- filter(dataset, Symbol %in% total.process.string)
  process.heatmap <- filter(TPM_HM, RefSeq %in% process.heatmap[[1]])
  process.input <- process.heatmap[, c(3:ncol(process.heatmap))]
  rownames(process.input) <- process.heatmap[,2]
  process.input <- t(process.input)
  
  Heatmap(process.input,
          cluster_columns = TRUE,
          column_names_side = "top", column_dend_side = "top",
          row_names_side = "left", row_dend_side = "left",
          row_names_gp = gpar(fontsize = 14), # fontface = "bold"
          row_dend_width = unit(2, "cm"),
          clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
          name = "z-score",
          # column_title = paste0(process$Description),
          # row_title = "Silencing conditions",
          column_title_gp = gpar(fontsize = 14, fontfamily = "Arial"), # fontface = "bold",
          row_title_rot = 0,
          col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
          rect_gp = gpar(col = "black"),
          row_title_gp = gpar(fontsize = 40), # fontface = "bold"
  heatmap_legend_param = list(
    title = "z-score", 
    title_gp = gpar(fontsize = 14), # fontface = "bold"
    title_position = "leftcenter-rot",
    direction = "vertical"
  ))
}

total.process.string <- c()
heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0006954", cutoff = 0.6)

saveRDS(total.process.string, file.path(proj.dir, data.dir, "inflammatory_response_genes.RDS"))

inflammatory_response_genes <- readRDS(file.path(proj.dir, data.dir, "inflammatory_response_genes.RDS"))
inflammatory_response_genes <- t(inflammatory_response_genes)
write.table(inflammatory_response_genes, file.path(proj.dir, data.dir, "inflammatory_response_genes.txt"),
            row.names = FALSE, quote = FALSE, sep = ", ")

total.process.string <- c()
heatmap.generator(dataset = DF, TPM_HM = TPM_HM, list.of.BP = BP.sim.cutoff, 
                  GOid = "GO:0040008", cutoff = 0.6)

saveRDS(total.process.string, file.path(proj.dir, data.dir, "growth_regulation_genes.RDS"))

growth_regulation_genes <- readRDS(file.path(proj.dir, data.dir, "growth_regulation_genes.RDS"))
growth_regulation_genes <- t(growth_regulation_genes)
write.table(growth_regulation_genes, file.path(proj.dir, data.dir, "growth_regulation_genes.txt"),
            row.names = FALSE, quote = FALSE, sep = ", ")

#### end of modification

flog.debug("Session Info")
sessionInfo()
