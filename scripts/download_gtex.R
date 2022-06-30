#'#################################################################################
#'#################################################################################
#' Download and preprocess GTEx data
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(recount)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)

genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)

options(timeout=100000)
## Prepare all GTEx
download_retry(url = "http://duffel.rail.bio/recount/v2/SRP012682/rse_gene.Rdata",
  destfile = 'data/GTEx/rse_gene_all.Rdata', mode = "wb")

load('data/GTEx/rse_gene_all.Rdata')


assay(rse_gene)[assay(rse_gene) > .Machine$integer.max] <- .Machine$integer.max
gtex <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = colData(rse_gene),
                              design = ~ 1)
rownames(gtex) <- gsub("\\.[0-9]*", "", rownames(gtex) , perl = TRUE)
gtex <- vst(gtex[genes, ], blind=FALSE)
saveHDF5SummarizedExperiment(gtex, "results/GTEx/", prefix = "vst_all_")

tissue <- gtex$smtsd 
write.table(tissue, file = "results/GTEx/individuals_labels.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)


## Prepare prostate data
download_retry(url = "http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_prostate.Rdata",
  destfile = 'data/GTEx/rse_gene_prostate.Rdata', mode = "wb")

load('data/GTEx/rse_gene_prostate.Rdata')

prostate <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = colData(rse_gene),
                              design = ~ 1)
rownames(prostate) <- gsub("\\.[0-9]*", "", rownames(prostate) , perl = TRUE)
prostate <- vst(prostate[genes, ], blind=FALSE)
saveHDF5SummarizedExperiment(prostate, "results/GTEx/", prefix = "vst_prostate_")



## Prepare testis data
download_retry(url = "http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_testis.Rdata",
  destfile = 'data/GTEx/rse_gene_testis.Rdata', mode = "wb")

load('data/GTEx/rse_gene_testis.Rdata')

testis <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = colData(rse_gene),
                              design = ~ 1)
rownames(testis) <- gsub("\\.[0-9]*", "", rownames(testis) , perl = TRUE)
testis <- vst(testis[genes, ], blind=FALSE)
saveHDF5SummarizedExperiment(testis, "results/GTEx/", prefix = "vst_testis_")
