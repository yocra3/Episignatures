#'#################################################################################
#'#################################################################################
#' Evaluate GTEx clustering
#'#################################################################################
#'#################################################################################


## Load libraries
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)

## Load gtex data
prostate <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_prostate_")
testis <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_testis_")

path.map <- read.table("results/preprocess/go_kegg_final_gene_map.tsv", header = TRUE)

prost.feat <- read.table("results/GTEx_prostate/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE).
testis.feat <- read.table("results/GTEx_testis/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)

paths <- read.table("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(prost.feat) <- colnames(testis.feat) <-paths.vec

## Compare PCs genes and pathways
pc_genes <- prcomp(t(cbind(assay(prostate), assay(testis))))

pc_genes_df <- data.frame(pc_genes$x, tissue = rep(c("Prostate", "Testis"), c(ncol(prostate), ncol(testis))))
ggplot(pc_genes_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  theme_bw()

pc_paths <- prcomp(rbind(prost.feat, testis.feat ))

pc_paths_df <- data.frame(pc_paths$x, tissue = rep(c("Prostate", "Testis"), c(nrow(prost.feat), nrow(testis.feat))))
ggplot(pc_paths_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  theme_bw()
