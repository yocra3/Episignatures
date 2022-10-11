#'#################################################################################
#'#################################################################################
#' Prepare PRAD datasets for use in classification
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(cowplot)
library(tidyverse)
library(SummarizedExperiment)
library(HDF5Array)
library(NetActivity)
library(sva)

## GSE46691
load("results/preprocess/GSE46691/GSE46691_scores.Rdata")
gse46691_scores$gleason <- ifelse(gse46691_scores$`gleason score:ch1` >= 8, "High", "Low")

## GSE141551
load("results/preprocess/GSE141551/GSE141551_scores.Rdata")
gse141551_scores$gleason <- ifelse(gse141551_scores$`gleason_group:ch1` == "8-10", "High", "Low")

## GSE21034 - array
load("results/preprocess/GSE21034/GSE21034_scores.Rdata")
gse21034_scores$gleason <- ifelse(gse21034_scores$`biopsy_gleason_grade:ch1` >= 8, "High", "Low")

## GSE70768 - array
load("results/preprocess/GSE70768/GSE70768_scores.Rdata")
gse70768_scores$gleason <- ifelse(sapply(strsplit(gse70768_scores$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low")

## GSE70769 - array
load("results/preprocess/GSE70769/GSE70769_scores.Rdata")
gse70769_scores$gleason <- ifelse(sapply(strsplit(gse70769_scores$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low")


## GSE183019 - RNAseq
load("results/preprocess/GSE183019/GSE183019_scores.Rdata")
gse183019_scores$gleason <- ifelse(gse183019_scores$`gleason score:ch1` %in% c("4+4", "4+5", "5+4", "5+5 T4"), "High", "Low")
gse183019_scores$age <- as.numeric(gse183019_scores$`age:ch1`)

mod_gse183019 <- model.matrix(~ gleason + age, colData(gse183019_scores))
mod0_gse183019 <- model.matrix(~ age, colData(gse183019_scores))

n.sv <- num.sv(assay(gse183019_scores), mod_gse183019, method="leek")
svobj <- sva(assay(gse183019_scores), mod_gse183019, mod0_gse183019, n.sv=n.sv)
clean <- residuals(lmFit(assay(gse183019_scores), svobj$sv), assay(gse183019_scores))
assay(gse183019_scores) <- clean



## PRAD - RNAseq
prad <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad_tumor")
prad$gleason <- ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")
assay(prad) <- data.matrix(assay(prad))

prad_prep <- prepareSummarizedExperiment(prad, "gtex_gokegg")
prad_scores <- computeGeneSetScores(prad_prep, "gtex_gokegg")


## GSE169038
gse169038_se <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")
assay(gse169038_se) <- data.matrix(assay(gse169038_se))
gse169038_se$primary <- gsub("primary gleason: ", "", gse169038_se$characteristics_ch1.1)
gse169038_se$secondary <- gsub("secondary gleason: ", "", gse169038_se$characteristics_ch1.2)
gse169038_se$primary <- as.numeric(ifelse(gse169038_se$primary == "--", 1, gse169038_se$primary))
gse169038_se$secondary <- as.numeric(ifelse(gse169038_se$secondary == "--", 1, gse169038_se$secondary))
gse169038_se$gleason_cat <- paste(gse169038_se$primary, gse169038_se$secondary, sep = "-")
gse169038_se$gleason <- ifelse(gse169038_se$primary == 5 |  gse169038_se$secondary == 5 | gse169038_se$gleason_cat == "4+4", "High", "Low")

gse169038_se_filt <- gse169038_se[, !(gse169038_se$primary == 1 |  gse169038_se$secondary == 1)]
gse169038_prep <- prepareSummarizedExperiment(gse169038_se_filt, "gtex_gokegg")
gse169038_scores <- computeGeneSetScores(gse169038_prep, "gtex_gokegg")

## GSE201284 - RNAseq
load("results/preprocess/GSE201284/GSE201284_scores.Rdata")
gse201284_scores$gleason <- ifelse(gse201284_scores$`patient gleason score:ch1` %in% c("3+5", "3+5 T4", "4+4", "4+5", "5+4", "5+5", "4+4 T5", "4+5 T3", "5+4 T3"), "High", "Low")

prad_meta <- read_table("results/PRAD_metanalysis/PRAD_metanalysis1.tbl")
prad_meta$p_adj <- p.adjust(prad_meta$`P-value`)
meta_sig <- subset(prad_meta, p_adj < 0.05)


### Combine training
removeColData <- function(SE){
  colData(SE) <- colData(SE)[, "gleason", drop = FALSE]
  SE
}
prad_train <- Reduce(cbind, lapply(list(gse70769_scores, gse21034_scores, gse183019_scores, gse141551_scores, gse201284_scores, gse169038_scores), removeColData))
prad_train$dataset <- rep(c("gse70769", "gse21034", "gse183019", "gse141551", "gse201284", "gse169038"),
                          sapply(list(gse70769_scores, gse21034_scores, gse183019_scores, gse141551_scores, gse201284_scores, gse169038_scores), ncol))

prad_train_pc <- prcomp(t(assay(prad_train)), rank. = 4)
prad_train_pc$x %>%
  data.frame() %>%
  mutate(dataset = prad_train$dataset, gleason = prad_train$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = dataset, shape = gleason)) +
  geom_point() +
  theme_bw()
summary(lm(prad_train_pc$x[, 1] ~ dataset + gleason, colData(prad_train)))
summary(lm(prad_train_pc$x[, 2] ~ dataset + gleason, colData(prad_train)))

prad_train_sel_pc <- prcomp(t(assay(prad_train[meta_sig$MarkerName, ])), rank. = 4)
prad_train_sel_pc$x %>%
  data.frame() %>%
  mutate(dataset = prad_train$dataset, gleason = prad_train$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = dataset, shape = gleason)) +
  geom_point() +
  theme_bw()
summary(lm(prad_train_sel_pc$x[, 1] ~ dataset + gleason, colData(prad_train)))
summary(lm(prad_train_sel_pc$x[, 2] ~ dataset + gleason, colData(prad_train)))



saveHDF5SummarizedExperiment(prad_train, "results/PRAD_metanalysis/", prefix = "prad_train_")
saveHDF5SummarizedExperiment(prad_train[meta_sig$MarkerName, ], "results/PRAD_metanalysis/", prefix = "prad_train_sel_")

write.table(colData(prad_train), file = "results/PRAD_metanalysis/individuals_labels.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
