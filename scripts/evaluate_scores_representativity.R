#'#################################################################################
#'#################################################################################
#' Explore whether using NetActivity can correct batch effect - FAILED
#'#################################################################################
#'#################################################################################

## Load libraries
library(SummarizedExperiment)
library(HDF5Array)
library(DESeq2)
library(NetActivity)
library(NetActivityData)
library(tidyverse)
library(cowplot)
library(sva)

## Load TCGA data
load("data/tcga_gexp_combat.Rdata")
data(gtex_gokegg)
## Prostate
prad <- gexp_tcga_combat[, gexp_tcga_combat$project_id == "TCGA-PRAD"]
load('data/GTEx/rse_gene_prostate.Rdata')
prostate <- rse_gene
rownames(prostate) <- gsub("\\.[0-9]*", "", rownames(prostate) , perl = TRUE)

load('data/GTEx/rse_gene_testis.Rdata')
testis <- rse_gene
rownames(testis) <- gsub("\\.[0-9]*", "", rownames(testis) , perl = TRUE)



com_genes <- intersect(rownames(prostate), rownames(prad))
se_prostate <- SummarizedExperiment(cbind(assay(prostate[com_genes, ]),
                                    assay(prad[com_genes, ])),
              colData = data.frame(batch = rep(c("GTEx", "TCGA"), c(ncol(prostate), ncol(prad))),
                      status = c(rep("Control", ncol(prostate)),
                                  ifelse(prad$sample_type == "Solid Tissue Normal", "Control", "Cancer"))),
              rowData = rowData(prad[com_genes, ]))

# Apply ComBat
comb_prostate <- ComBat_seq(assay(se_prostate), batch = se_prostate$batch , group = se_prostate$status, full_mod=TRUE)
deseq_prostate <- DESeqDataSetFromMatrix(
                    countData = assay(se_prostate),
                    colData = colData(se_prostate),
                    design = ~ batch + status)
vst_prostate <-  vst(deseq_prostate, blind=FALSE)

pc_prost_vst <- prcomp(t(assay(vst_prostate)), rank. = 4)

pc_prost_vst_plot <- pc_prost_vst$x %>%
  as_tibble() %>%
  mutate(Status = se_prostate$status,
         Batch = se_prostate$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Status, shape = Batch)) +
  geom_point() +
  theme_bw() + ggtitle("Vst prostate")

preproc_prostate <- prepareSummarizedExperiment(vst_prostate, "gtex_gokegg")
pc_prost_preproc <- prcomp(t(assay(preproc_prostate[colnames(gtex_gokegg), ])), rank. = 4)
pc_prost_preproc_plot <- pc_prost_preproc$x %>%
  as_tibble() %>%
  mutate(Status = se_prostate$status,
         Batch = se_prostate$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Status, shape = Batch)) +
  geom_point() +
  theme_bw() +
  ggtitle("Preproc Prostate")

scores_prostate <- computeGeneSetScores(preproc_prostate, "gtex_gokegg")
pc_scores_preproc <- prcomp(t(assay(scores_prostate)), rank. = 4)
pc_prost_score_plot <- pc_scores_preproc$x %>%
  as_tibble() %>%
  mutate(Status = se_prostate$status,
         Batch = se_prostate$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Status, shape = Batch)) +
  geom_point() +
  theme_bw() +
  ggtitle("Scores Prostate")

plot_grid(pc_prost_vst_plot, pc_prost_preproc_plot, pc_prost_score_plot, nrow = 1)


## Combine two GTEx tissues
se_mix <- SummarizedExperiment(cbind(assay(testis[com_genes, ]),
                                    assay(prostate[com_genes, ]),
                                    assay(prad[com_genes, prad$sample_type == "Solid Tissue Normal"])),
              colData = data.frame(batch = rep(c("GTEx", "TCGA"), c(ncol(prostate) + ncol(testis), sum(prad$sample_type == "Solid Tissue Normal"))),
                      tissue = rep(c("Testis", "Prostate"), c(ncol(testis), ncol(prostate) + sum(prad$sample_type == "Solid Tissue Normal")))),
              rowData = rowData(prad[com_genes, ]))
#
deseq_mix <- DESeqDataSetFromMatrix(
                    countData = assay(se_mix),
                    colData = colData(se_mix),
                    design = ~ batch + tissue)
vst_mix <-  vst(deseq_mix, blind=FALSE)

pc_mix_vst <- prcomp(t(assay(vst_mix)), rank. = 4)

pc_mix_vst_plot <- pc_mix_vst$x %>%
  as_tibble() %>%
  mutate(Tissue = deseq_mix$tissue,
         Batch = deseq_mix$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Tissue, shape = Batch)) +
  geom_point() +
  theme_bw() + ggtitle("Vst Mix")

#

preproc_mix <- prepareSummarizedExperiment(vst_mix, "gtex_gokegg")
pc_mix_preproc <- prcomp(t(assay(preproc_mix[colnames(gtex_gokegg), ])), rank. = 4)
pc_mix_preproc_plot <- pc_mix_preproc$x %>%
  as_tibble() %>%
  mutate(Tissue = deseq_mix$tissue,
         Batch = deseq_mix$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Tissue, shape = Batch)) +
  geom_point() +
  theme_bw() +
  ggtitle("Preproc mix")

scores_mix <- computeGeneSetScores(preproc_mix, "gtex_gokegg")
pc_scores_preproc <- prcomp(t(assay(scores_mix)), rank. = 4)
pc_mix_score_plot <- pc_scores_preproc$x %>%
  as_tibble() %>%
  mutate(Tissue = deseq_mix$tissue,
         Batch = deseq_mix$batch) %>%
  ggplot(aes(x = PC1, y = PC2, color = Tissue, shape = Batch)) +
  geom_point() +
  theme_bw() +
  ggtitle("Scores mix")

plot_grid(pc_mix_vst_plot, pc_mix_preproc_plot, pc_mix_score_plot, nrow = 1)
