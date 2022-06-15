#'#################################################################################
#'#################################################################################
#' Process HELIX gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(hta20transcriptcluster.db)
library(SummarizedExperiment)
library(HDF5Array)
library(tidyverse)

## Cobert IDs
load("data/HELIX/genexpr.Rdata")

mapping <- mapIds(
  hta20transcriptcluster.db,
  keys = rownames(genexpr),
  column = 'ENSEMBL',
  keytype = 'PROBEID')

se <- SummarizedExperiment(exprs(genexpr), colData = pData(genexpr), rowData = fData(genexpr))
rowData(se)$gene <- mapping
save(se, file = "results/HELIX/allGenes.se.RData")

## Subset to genes present in TCGA
genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)

## For each gene, select TC with more probes and higher Call rate
sel_TCs <- rowData(se) %>%
  data.frame() %>%
  group_by(gene) %>%
  slice_max(total_probes) %>%
  slice_max(CallRate) %>%
  slice_head(n = 1) %>%
  filter(gene %in% genes)

se.filt <- se[sel_TCs$transcript_cluster_id , ]
rownames(se.filt) <- rowData(se.filt)$gene

out_probes <- setdiff(genes, rownames(se.filt))
out <- matrix(0, ncol(se.filt),
              nrow = length(out_probes), ncol = ncol(se.filt),
              dimnames = list(out_probes, colnames(se.filt)))

new_assay <- rbind(assay(se.filt), out)
se.tcga_genes <- SummarizedExperiment(new_assay, colData = colData(se),
  rowData = data.frame(transcript_id = c(rowData(se.filt)$transcript_cluster_id, rep("Not present", length(out_probes)))))
se.tcga_genes <- se.tcga_genes[genes, ]
saveHDF5SummarizedExperiment(se.tcga_genes, "results/HELIX/", prefix = "network_genes")
