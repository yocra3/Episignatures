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

genes <- as.character(read.table("data/gene_names_GEX.csv")$V1)
load('data/GTEx/rse_gene_all.Rdata')


assay(rse_gene)[assay(rse_gene) > .Machine$integer.max] <- .Machine$integer.max
gtex <- DESeqDataSetFromMatrix(countData = assay(rse_gene),
                              colData = colData(rse_gene),
                              design = ~ smtsd)
rownames(gtex) <- gsub("\\.[0-9]*", "", rownames(gtex) , perl = TRUE)

missing <-  genes[which(!genes %in% rownames(gtex))]
write.table(missing, file = "results/Kike/missing_genes.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

gtex_filt <- gtex[genes[!genes %in% missing], ]
gtex_vst <- vst(gtex_filt, blind=FALSE)
saveHDF5SummarizedExperiment(gtex_vst, "results/Kike/", prefix = "vst_all_")

tissue <- gtex_vst$smtsd
write.table(tissue, file = "results/Kike/individuals_labels.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

#
sel_genes <- rownames(gtex_vst)
write.table(sel_genes, file = "results/Kike/input_genes.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
