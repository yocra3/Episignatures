#'#################################################################################
#'#################################################################################
#' Download and preprocess GSE54460 data
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(recount)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(GEOquery)

genes <- as.character(read.table("./results/GTEx_coding/input_genes.txt")$V1)

options(timeout=100000)
## Prepare all GTEx
download_retry(url = "http://duffel.rail.bio/recount/v2/SRP036848/rse_gene.Rdata",
  destfile = 'data/GSE54460/rse_gene_all.Rdata', mode = "wb")

load('data/GSE54460/rse_gene_all.Rdata')

## Download GEO
geo <- getGEO("GSE54460")

assay(rse_gene)[assay(rse_gene) > .Machine$integer.max] <- .Machine$integer.max
colnames(rse_gene) <- rse_gene$geo_accession
gse54460 <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = pData(geo[[1]][, colnames(rse_gene)]),
                              design = ~ 1)
rownames(gse54460) <- gsub("\\.[0-9]*", "", rownames(gse54460) , perl = TRUE)
gse54460 <- vst(gse54460[genes, ], blind=FALSE)
saveHDF5SummarizedExperiment(gse54460, "results/GSE54460/", prefix = "vst_all_")
