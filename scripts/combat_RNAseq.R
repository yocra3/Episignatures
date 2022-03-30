#! /usr/local/bin/Rscript

#'#################################################################################
#'#################################################################################
#' Process TCGA gene expression for deep learning
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(SummarizedExperiment)
library(tidyverse)
library(HDF5Array)
library(sva)

## Load data
load("data/tcga_gexp.Rdata")
gexp_fold <- "results/TCGA_gexp/"

### Correct sequencing center with ComBat

get_IDS <- function (data)   {
IDs <- strsplit(c(colnames(data)), "-")
IDs <- plyr::ldply(IDs, rbind)
colnames(IDs) <- c("project", "tss", "participant", "sample", "portion", "plate", "center")
cols <- c("project", "tss", "participant")
IDs$patient <- apply(IDs[, cols], 1, paste, collapse = "-")
barcode <- colnames(data)
IDs <- cbind(IDs, barcode)
condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
condition <- gsub("01+[[:alpha:]]", "cancer", condition)
IDs$condition <- condition
IDs$myorder <- 1:nrow(IDs)
return(IDs)
}
phenos <- get_IDS(gexp_tcga)

gexp_tcga$sample_type2 <- ifelse(gexp_tcga$sample_type == "Solid Tissue Normal", "normal", "tumor")

adj.mod <- model.matrix(~ project_id, colData(gexp_tcga))[, -14]## Avoid confounding
adj_counts <- ComBat_seq(assay(gexp_tcga), batch = phenos$center , group = gexp_tcga$sample_type2, covar_mod = adj.mod, full_mod=TRUE)
save(adj_counts, file = "data/combat_counts.Rdata")
