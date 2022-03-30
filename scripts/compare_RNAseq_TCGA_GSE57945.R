#'#################################################################################
#'#################################################################################
#' Process TCGA gene expression for deep learning
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(sva)
library(rtracklayer)

## Load data
load("data/tcga_gexp_combat.Rdata")
vst.tcga <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat/", prefix = "vsd_norm")
load( "results/SRP042228/RSE_phenotypes.Rdata")
vst.gse <- loadHDF5SummarizedExperiment( "results/SRP042228/", prefix = "vsd_norm")
rownames(rse_gene) <- gsub("\\.[0-9]*", "", rownames(rse_gene), perl = TRUE)

## Load gencode annotation
annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
annot <- filter(annot, type == "gene") %>%
  as_tibble()
annot_cod <- subset(annot, gene_type == "protein_coding")   %>%
   mutate(ensembl = gsub("\\.[0-9]*_[0-9]*", "", gene_id , perl = TRUE))

### Correct sequencing center with ComBat
phenos <- get_IDs(gexp_tcga_combat)

### Compute medians for original vst
tcga_medians <- rowMedians(data.matrix(assay(vst.tcga)))
ctrl_medians <- rowMedians(data.matrix(assay(vst.tcga[, vst.tcga$sample_type == "Solid Tissue Normal"])))
names(ctrl_medians) <- names(tcga_medians) <- rownames(vst.tcga)

gse_medians <- rowMedians(data.matrix(assay(vst.gse)))
names(gse_medians) <- rownames(vst.gse)

com_genes <- intersect(names(ctrl_medians) , names(gse_medians) )

png("figures/rnaseq_comp/vst_ori_medians.png", width = 1500)
par(mfrow = c(1, 3))
plot(tcga_medians[com_genes], gse_medians[com_genes], main = "All TCGA vs GSE")
plot(ctrl_medians[com_genes], gse_medians[com_genes], main = "Control TCGA vs GSE")
plot(tcga_medians[com_genes], ctrl_medians[com_genes], main = "All TCGA vs Control TCGA")
dev.off()

summary(lm(tcga_medians[com_genes] ~ gse_medians[com_genes]))
summary(lm(ctrl_medians[com_genes] ~ gse_medians[com_genes]))
summary(lm(tcga_medians[com_genes] ~ ctrl_medians[com_genes]))

### Recompute VST including only coding genes
gexp_tcga_coding <- gexp_tcga_combat[rownames(gexp_tcga_combat) %in% annot_cod$ensembl, ]

adj_counts2 <- assay(gexp_tcga_coding)
mode(adj_counts2) <- "integer"

dds_tcga <- DESeqDataSetFromMatrix(countData = assay(gexp_tcga_coding[rowSums(is.na(adj_counts2)) == 0, ]),
                              colData = colData(gexp_tcga_coding),
                              design = ~ project_id)
vst_tcga_cod <- vst(dds_tcga, blind=FALSE)

rse_gene_coding <- rse_gene[rownames(rse_gene) %in% annot_cod$ensembl, ]
dds_gse <- DESeqDataSetFromMatrix(countData = assay(rse_gene_coding),
                              colData = colData(rse_gene_coding),
                              design = ~ 1)
vst_gse_cod <- vst(dds_gse, blind=FALSE)


### Compute medians for vst of coding genes
tcga_cod_medians <- rowMedians(data.matrix(assay(vst_tcga_cod)))
ctrl_cod_medians <- rowMedians(data.matrix(assay(vst_tcga_cod[, vst_tcga_cod$sample_type == "Solid Tissue Normal"])))
names(ctrl_cod_medians) <- names(tcga_cod_medians) <- rownames(vst_tcga_cod)

gse_cod_medians <- rowMedians(data.matrix(assay(vst_gse_cod)))
names(gse_cod_medians) <- rownames(vst_gse_cod)

com_genes_cod <- intersect(names(ctrl_cod_medians) , names(tcga_cod_medians) )

png("figures/rnaseq_comp/vst_coding_medians.png", width = 1500)
par(mfrow = c(1, 3))
plot(tcga_cod_medians[com_genes_cod], gse_cod_medians[com_genes_cod], main = "All TCGA vs GSE")
plot(ctrl_cod_medians[com_genes_cod], gse_cod_medians[com_genes_cod], main = "Control TCGA vs GSE")
plot(tcga_cod_medians[com_genes_cod], ctrl_cod_medians[com_genes_cod], main = "All TCGA vs Control TCGA")
dev.off()

summary(lm(tcga_cod_medians[com_genes_cod] ~ gse_cod_medians[com_genes_cod]))
summary(lm(ctrl_cod_medians[com_genes_cod] ~ gse_cod_medians[com_genes_cod]))
summary(lm(tcga_cod_medians[com_genes_cod] ~ ctrl_cod_medians[com_genes_cod]))


### Apply ComBat seq
merge_rse <- SummarizedExperiment(assays = cbind(assay(rse_gene[com_genes, ]), assay(gexp_tcga_combat[com_genes, gexp_tcga_combat$sample_type == "Solid Tissue Normal"])),
                                  rowRanges = rowRanges(rse_gene[com_genes, ]))
merge_rse$project <- ifelse(grepl("TCGA", colnames(merge_rse)), "TCGA", "GSE")
merge_combat <- ComBat_seq(assay(merge_rse), batch = merge_rse$project)

merge_combat2 <- merge_combat
mode(merge_combat2) <- "integer"
merge_combat2[is.na(merge_combat2)] <- max(merge_combat2, na.rm = T)

vst_ctrl_combat <- DESeqDataSetFromMatrix(countData = merge_combat2[,merge_rse$project == "TCGA"],
                                          colData = colData(merge_rse[,merge_rse$project == "TCGA"]),
                                          design = ~ 1) %>%
                                          vst(blind=FALSE)
vst_gse_combat <- DESeqDataSetFromMatrix(countData = merge_combat2[,merge_rse$project == "GSE"],
                                          colData = colData(merge_rse[,merge_rse$project == "GSE"]),
                                          design = ~ 1) %>%
                                          vst(blind=FALSE)

### Compute medians for vst combat
ctrl_combat_medians <- rowMedians(data.matrix(assay(vst_ctrl_combat)))
gse_combat_medians <- rowMedians(data.matrix(assay(vst_gse_combat)))
names(ctrl_combat_medians) <- names(gse_combat_medians) <- rownames(vst_ctrl_combat)

png("figures/rnaseq_comp/vst_combat_medians.png", width = 1500)
par(mfrow = c(1, 3))
plot(ctrl_combat_medians, gse_combat_medians, main = "Control TCGA ComBat vs GSE ComBat",
  col = ifelse(names(ctrl_combat_medians) %in% annot_cod$ensembl, "red", "black"))
plot(ctrl_medians[com_genes], ctrl_combat_medians[com_genes], main = "Original vs Combat - Control TCGA",
  col = ifelse(com_genes %in% annot_cod$ensembl, "red", "black"))
plot(ctrl_medians[com_genes], gse_combat_medians[com_genes], main = "Control TCGA original vs GSE ComBat",
  col = ifelse(com_genes %in% annot_cod$ensembl, "red", "black"))
dev.off()

summary(lm(ctrl_combat_medians ~ gse_combat_medians)
summary(lm(ctrl_medians[com_genes] ~ ctrl_combat_medians[com_genes]))
summary(lm(tcga_medians[com_genes] ~ gse_combat_medians[com_genes]))

summary(lm(ctrl_combat_medians ~ gse_combat_medians, subset = names(ctrl_combat_medians) %in% annot_cod$ensembl))
summary(lm(ctrl_medians[com_genes] ~ ctrl_combat_medians[com_genes], subset = com_genes %in% annot_cod$ensembl))
summary(lm(tcga_medians[com_genes] ~ gse_combat_medians[com_genes], subset = com_genes %in% annot_cod$ensembl))


saveHDF5SummarizedExperiment(vst_ctrl_combat, "figures/rnaseq_comp/", prefix = "ctrl_combat_")
saveHDF5SummarizedExperiment(vst_gse_combat, "figures/rnaseq_comp/", prefix = "gse_combat_")


pc_merge <- prcomp(t(cbind(assay(vst_ctrl_combat), assay(vst_gse_combat))))
png("figures/rnaseq_comp/vst_combat_pca.png")
plot(pc_merge$x, col = rep(c("black", "red"), c(ncol(vst_ctrl_combat), ncol(vst_gse_combat))))
dev.off()

### Apply ComBat to vst
merge_vst_rse <- SummarizedExperiment(assays = cbind(assay(vst.gse[com_genes, ]), assay(vst.tcga[com_genes, vst.tcga$sample_type == "Solid Tissue Normal"])),
                                  rowRanges = rowRanges(vst.gse[com_genes, ]))
merge_vst_rse$project <- ifelse(grepl("TCGA", colnames(merge_vst_rse)), "TCGA", "GSE")
merge_vst_combat <- ComBat(assay(merge_vst_rse), batch = merge_vst_rse$project)

vst2_ctrl_combat <- merge_vst_combat[, merge_vst_rse$project == "TCGA"]
vst2_gse_combat <- merge_vst_combat[, merge_vst_rse$project == "GSE"]


### Compute medians for vst combat
ctrl_combat2_medians <- rowMedians(vst2_ctrl_combat)
gse_combat2_medians <- rowMedians(vst2_gse_combat)
names(ctrl_combat2_medians) <- names(gse_combat2_medians) <- rownames(vst2_ctrl_combat)

png("figures/rnaseq_comp/vst_combat2_medians.png", width = 1500)
par(mfrow = c(1, 3))
plot(ctrl_combat2_medians, gse_combat2_medians, main = "Control TCGA ComBat vs GSE ComBat",
  col = ifelse(names(ctrl_combat_medians) %in% annot_cod$ensembl, "red", "black"))
plot(ctrl_medians[com_genes], ctrl_combat2_medians[com_genes], main = "Original vs Combat - Control TCGA",
  col = ifelse(com_genes %in% annot_cod$ensembl, "red", "black"))
plot(ctrl_medians[com_genes], gse_combat2_medians[com_genes], main = "Control TCGA original vs GSE ComBat",
  col = ifelse(com_genes %in% annot_cod$ensembl, "red", "black"))
dev.off()

summary(lm(ctrl_combat2_medians ~ gse_combat2_medians))
summary(lm(ctrl_medians[com_genes] ~ ctrl_combat2_medians[com_genes]))
summary(lm(tcga_medians[com_genes] ~ gse_combat2_medians[com_genes]))

pc_merge2 <- prcomp(t(merge_vst_combat))
png("figures/rnaseq_comp/vst_combat2_pca.png")
plot(pc_merge2$x, col = rep(c("black", "red"), col = factor(merge_vst_rse$project)))
dev.off()

vst_ctrl_combat2 <- SummarizedExperiment(assays = vst2_ctrl_combat,
                                  rowRanges = rowRanges(merge_vst_rse))
vst_gse_combat2 <- SummarizedExperiment(assays = vst2_gse_combat,
                                  rowRanges = rowRanges(merge_vst_rse))
saveHDF5SummarizedExperiment(vst_ctrl_combat2, "figures/rnaseq_comp/", prefix = "ctrl_combat2_")
saveHDF5SummarizedExperiment(vst_gse_combat2, "figures/rnaseq_comp/", prefix = "gse_combat2_")
