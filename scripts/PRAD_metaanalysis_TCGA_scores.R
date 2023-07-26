#'#################################################################################
#'#################################################################################
#' Run a meta-analysis in PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(cowplot)
library(tidyverse)
library(SummarizedExperiment)
library(HDF5Array)
library(NetActivity)
library(cowplot)

createMETAL_tab <- function(fit){

  SE <- sqrt(fit$s2.post) * fit$stdev.unscaled %>% data.frame()
  res <- topTable(fit, coef = 2, n = Inf)
  res$SE <- SE[rownames(res), 2]
  res$Marker <- rownames(res)
  res$TestedAllele <- "C"
  res$OtherAllele <- "T"
  res$N <- nrow(fit$qr$qr)
  res
}
dir.create("results/PRAD_metanalysis_TCGA/")

## GSE46691
load("results/preprocess/GSE46691/GSE46691_scores_tcga.Rdata")

gse46691_scores_tcga$gleason <- factor(ifelse(gse46691_scores_tcga$`gleason score:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse46691_scores_tcga$metastasis <- gse46691_scores_tcga$`metastatic event:ch1`

mod_gse46691 <- model.matrix(~ gleason + metastasis, colData(gse46691_scores_tcga))
lm.gse46691 <- lmFit(assay(gse46691_scores_tcga), mod_gse46691) %>% eBayes()
tab.gse46691 <- topTable(lm.gse46691, coef = 2, n = Inf)

metal.gse46691 <- createMETAL_tab(lm.gse46691)
write.table(metal.gse46691, file = "results/PRAD_metanalysis_TCGA/GSE46691_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

gse46691_pc <- prcomp(t(assay(gse46691_scores_tcga)), rank. = 4)
gse46691_pc_plot <- gse46691_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse46691_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE46691") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE141551
load("results/preprocess/GSE141551/GSE141551_scores_tcga.Rdata")
gse141551_scores_tcga$gleason <- factor(ifelse(gse141551_scores_tcga$`gleason_group:ch1` == "8-10", "High", "Low"), levels = c("Low", "High"))
gse141551_scores_tcga$race <- gse141551_scores_tcga$`race:ch1`
gse141551_scores_tcga$age <- gse141551_scores_tcga$`age_group:ch1`
gse141551_scores_tcga$stage <- gse141551_scores_tcga$`stage_disease:ch1`

mod_gse141551 <- model.matrix(~ gleason + race + age + stage, colData(gse141551_scores_tcga))
lm.gse141551 <- lmFit(assay(gse141551_scores_tcga), mod_gse141551) %>% eBayes()
tab.gse141551 <- topTable(lm.gse141551, coef = 2, n = Inf)

metal.gse141551 <- createMETAL_tab(lm.gse141551)
write.table(metal.gse141551, file = "results/PRAD_metanalysis_TCGA/GSE141551_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse141551_pc <- prcomp(t(assay(gse141551_scores_tcga)), rank. = 4)
gse141551_pc_plot <- gse141551_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse141551_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE141551") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE21034 - array
load("results/preprocess/GSE21034/GSE21034_scores_tcga.Rdata")
gse21034_scores_tcga$gleason <- factor(ifelse(gse21034_scores_tcga$`biopsy_gleason_grade:ch1` >= 8, "High", "Low"), levels = c("Low", "High"))
gse21034_scores_tcga$metastasis <- gse21034_scores_tcga$`tumor type:ch1`

mod_gse21034 <- model.matrix(~ gleason + metastasis, colData(gse21034_scores_tcga))
lm.gse21034 <- lmFit(assay(gse21034_scores_tcga), mod_gse21034) %>% eBayes()
tab.gse21034 <- topTable(lm.gse21034, coef = 2, n = Inf)

metal.gse21034 <- createMETAL_tab(lm.gse21034)
write.table(metal.gse21034, file = "results/PRAD_metanalysis_TCGA/GSE21034_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse21034_pc <- prcomp(t(assay(gse21034_scores_tcga)), rank. = 4)
gse21034_pc_plot <- gse21034_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse21034_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE21034") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE70768 - array
load("results/preprocess/GSE70768/GSE70768_scores_tcga.Rdata")
gse70768_scores_tcga$gleason <- factor(ifelse(sapply(strsplit(gse70768_scores_tcga$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))
gse70768_scores_tcga$age <- as.numeric(gse70768_scores_tcga$`age at diag:ch1`)
gse70768_scores_tcga$tum_prop <- as.numeric(gsub("%", "", gse70768_scores_tcga$`tumour %:ch1`))

mod_gse70768 <- model.matrix(~ gleason + age + tum_prop, colData(gse70768_scores_tcga))
lm.gse70768 <- lmFit(assay(gse70768_scores_tcga)[, rownames(mod_gse70768)], mod_gse70768) %>% eBayes()
tab.gse70768 <- topTable(lm.gse70768, coef = 2, n = Inf)

metal.gse70768 <- createMETAL_tab(lm.gse70768)
write.table(metal.gse70768, file = "results/PRAD_metanalysis_TCGA/GSE70768_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse70768_pc <- prcomp(t(assay(gse70768_scores_tcga[, colSums(is.na(assay(gse70768_scores_tcga))) == 0])), rank. = 4)
gse70768_pc_plot <- gse70768_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse70768_scores_tcga[, colSums(is.na(assay(gse70768_scores_tcga))) == 0]$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE70768") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## GSE70769 - array
load("results/preprocess/GSE70769/GSE70769_scores_tcga.Rdata")
gse70769_scores_tcga$gleason <- factor(ifelse(sapply(strsplit(gse70769_scores_tcga$`tumour gleason:ch1`, "="), `[`, 1) %>% as.numeric() >= 8, "High", "Low"), levels = c("Low", "High"))

mod_gse70769 <- model.matrix(~ gleason, colData(gse70769_scores_tcga))
lm.gse70769 <- lmFit(assay(gse70769_scores_tcga), mod_gse70769) %>% eBayes()
tab.gse70769 <- topTable(lm.gse70769, coef = 2, n = Inf)

metal.gse70769 <- createMETAL_tab(lm.gse70769)
write.table(metal.gse70769, file = "results/PRAD_metanalysis_TCGA/GSE70769_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse70769_pc <- prcomp(t(assay(gse70769_scores_tcga)), rank. = 4)
gse70769_pc_plot <- gse70769_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse70769_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE70769") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



## GSE183019 - RNAseq
load("results/preprocess/GSE183019/GSE183019_scores_tcga.Rdata")
gse183019_scores_tcga$gleason <- factor(ifelse(gse183019_scores_tcga$`gleason score:ch1` %in% c("4+4", "4+5", "5+4", "5+5 T4"), "High", "Low"), levels = c("Low", "High"))
gse183019_scores_tcga$age <- as.numeric(gse183019_scores_tcga$`age:ch1`)


mod_gse183019 <- model.matrix(~ gleason + age, colData(gse183019_scores_tcga))
lm.gse183019 <- lmFit(assay(gse183019_scores_tcga), mod_gse183019) %>% eBayes()
tab.gse183019 <- topTable(lm.gse183019, coef = 2, n = Inf)

metal.gse183019 <- createMETAL_tab(lm.gse183019)
write.table(metal.gse183019, file = "results/PRAD_metanalysis_TCGA/GSE183019_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

n.sv <- num.sv(assay(gse183019_scores_tcga), mod_gse183019, method="leek")
mod0_gse183019 <- model.matrix(~ age, colData(gse183019_scores_tcga))
svobj <- sva(assay(gse183019_scores_tcga), mod_gse183019, mod0_gse183019, n.sv=n.sv)
modSv <- cbind(mod_gse183019, svobj$sv)
lm.gse183019_sva <- lmFit(assay(gse183019_scores_tcga), modSv) %>% eBayes()
tab.gse183019_sva <- topTable(lm.gse183019_sva, coef = 2, n = Inf)


metal.gse183019_sva <- createMETAL_tab(lm.gse183019_sva)
write.table(metal.gse183019_sva, file = "results/PRAD_metanalysis_TCGA/GSE183019_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse183019_pc <- prcomp(t(assay(gse183019_scores_tcga)), rank. = 4)
gse183019_pc_plot <- gse183019_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse183019_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE183019") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## PRAD - RNAseq
load("results/TCGA_PRAD/vst_SE.Rdata")

prad_prep <- prepareSummarizedExperiment(vst.prad, "tcga_gokegg")
prad_scores <- computeGeneSetScores(prad_prep, "tcga_gokegg")

mod_prad <- model.matrix(~ gleason + paper_Subtype + age_at_index + race, colData(vst.prad))
lm.prad <- lmFit(assay(prad_scores), mod_prad) %>% eBayes()

metal.prad <- createMETAL_tab(lm.prad)
write.table(metal.prad, file = "results/PRAD_metanalysis_TCGA/PRAD_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
prad_pc <- prcomp(t(assay(prad_scores)), rank. = 4)
prad_pc_plot <- prad_pc$x %>%
  data.frame() %>%
  mutate(Gleason = prad_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("PRAD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE169038
gse169038_se <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")
assay(gse169038_se) <- data.matrix(assay(gse169038_se))
gse169038_se$race <- ifelse(grepl( "White", gse169038_se$characteristics_ch1.4), "EUR", "AFR")
gse169038_se$decipher <- factor(gsub("decipher risk group: ", "", gse169038_se$characteristics_ch1.3), levels = c("Lower", "Average", "Higher"))
gse169038_se$primary <- gsub("primary gleason: ", "", gse169038_se$characteristics_ch1.1)
gse169038_se$secondary <- gsub("secondary gleason: ", "", gse169038_se$characteristics_ch1.2)
gse169038_se$primary <- as.numeric(ifelse(gse169038_se$primary == "--", 1, gse169038_se$primary))
gse169038_se$secondary <- as.numeric(ifelse(gse169038_se$secondary == "--", 1, gse169038_se$secondary))
gse169038_se$gleason_cat <- paste(gse169038_se$primary, gse169038_se$secondary, sep = "-")
gse169038_se$gleason <- factor(ifelse(gse169038_se$primary == 5 |  gse169038_se$secondary == 5 | gse169038_se$gleason_cat == "4+4", "High", "Low"), levels = c("Low", "High"))

gse169038_se_filt <- gse169038_se[, !(gse169038_se$primary == 1 |  gse169038_se$secondary == 1)]
gse169038_prep <- prepareSummarizedExperiment(gse169038_se_filt, "tcga_gokegg")
gse169038_scores <- computeGeneSetScores(gse169038_prep, "tcga_gokegg")

mod_gse169038 <- model.matrix(~ gleason + race + decipher, colData(gse169038_scores))
lm.gse169038 <- lmFit(assay(gse169038_scores), mod_gse169038) %>% eBayes()

metal.gse169038 <- createMETAL_tab(lm.gse169038)
write.table(metal.gse169038, file = "results/PRAD_metanalysis_TCGA/GSE169038_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse169038_pc <- prcomp(t(assay(gse169038_scores)), rank. = 4)
gse169038_pc_plot <- gse169038_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse169038_scores$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE169038") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## GSE201284 - RNAseq
load("results/preprocess/GSE201284/GSE201284_scores_tcga.Rdata")
gse201284_scores_tcga$gleason <- factor(ifelse(gse201284_scores_tcga$`patient gleason score:ch1` %in% c("3+5", "3+5 T4", "4+4", "4+5", "5+4", "5+5", "4+4 T5", "4+5 T3", "5+4 T3"), "High", "Low"), levels = c("Low", "High"))
gse201284_scores_tcga$age <- as.numeric(gse201284_scores_tcga$`age:ch1`)
gse201284_scores_tcga$tumor_prop <- as.numeric(gse201284_scores_tcga$`% tumor:ch1`)


mod_gse201284 <- model.matrix(~ gleason + age + tumor_prop, colData(gse201284_scores_tcga))
lm.gse201284 <- lmFit(assay(gse201284_scores_tcga)[, rownames(mod_gse201284)], mod_gse201284) %>% eBayes()
tab.gse201284 <- topTable(lm.gse201284, coef = 2, n = Inf)

metal.gse201284 <- createMETAL_tab(lm.gse201284)
write.table(metal.gse201284, file = "results/PRAD_metanalysis_TCGA/GSE201284_metal.txt",
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### PC
gse201284_pc <- prcomp(t(assay(gse201284_scores_tcga)), rank. = 4)
gse201284_pc_plot <- gse201284_pc$x %>%
  data.frame() %>%
  mutate(Gleason = gse201284_scores_tcga$gleason) %>%
  ggplot(aes(x = PC1, y = PC2, color = Gleason)) +
  geom_point() +
  theme_bw() +
  ggtitle("GSE201284") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## PCs of scores
png("figures/PRAD_metanalysis_TCGA_pcs.png", width = 1000, height = 1000)
plot_grid(prad_pc_plot, gse201284_pc_plot, gse183019_pc_plot,
          gse70769_pc_plot, gse70768_pc_plot, gse141551_pc_plot,
          gse46691_pc_plot, gse21034_pc_plot, gse169038_pc_plot,
        nrow = 3, labels = LETTERS[1:9])
dev.off()

## Evaluation
prad_meta <- read_table("results/PRAD_metanalysis/PRAD_metanalysis1.tbl")
prad_meta$p_adj <- p.adjust(prad_meta$`P-value`)

prad_meta_tcga <- read_table("results/PRAD_metanalysis_TCGA/PRAD_metanalysis1.tbl")
prad_meta_tcga$p_adj <- p.adjust(prad_meta_tcga$`P-value`)

dir_df <- strsplit(prad_meta_tcga$Direction, "") %>% Reduce(f = rbind) %>%
  data.frame()
colnames(dir_df) <- paste0("Dir_", c("GSE183019", "GSE70769", "GSE70768", "GSE46691",
  "GSE141551", "GSE21034", "PRAD", "GSE169038", "GSE201284"))

prad_meta_tcga <- cbind(prad_meta_tcga, dir_df) %>% as_tibble()
prad_meta_tcga <- mutate(prad_meta_tcga,
    pos = rowSums(across(starts_with("Dir_")) == "+"),
    concor = ifelse(pos < 5, 9 - pos, pos),
    com_dir = ifelse(sign(Effect) == -1, "-", "+"))



prad_meta_txt <- select(prad_meta_tcga, MarkerName, Direction, Effect, StdErr,  `P-value`, p_adj)
library(NetActivityData)
data(tcga_gokegg_annot)
colnames(prad_meta_txt) <- c("GeneSet", "Direction", "Effect", "StdErr", "pvalue", "adj_pvalue")
prad_meta_txt <- right_join(select(tcga_gokegg_annot, GeneSet, Term), prad_meta_txt, by = "GeneSet") %>%
  as_tibble() %>% arrange(pvalue)

write.table(prad_meta_txt, file = "results/PRAD_metanalysis_TCGA/PRAD_metanalysis.tsv",
          sep = "\t", quote = FALSE, row.names = FALSE)


ggplot(prad_meta, aes(x = factor(concor), y = -log10(`P-value`))) +
 geom_boxplot() +
 theme_bw()


meta_sig <- subset(prad_meta, p_adj < 0.05)
meta_sig_tcga <- subset(prad_meta_tcga, p_adj < 0.05)

meta_sig_comb <- full_join(select(prad_meta, MarkerName, Effect, `P-value`, p_adj),
                          select(prad_meta_tcga, MarkerName, Effect, `P-value`, p_adj),
                          suffix = c(".gtex", ".tcga"), by = "MarkerName") %>%
                  filter(p_adj.gtex < 0.05 | p_adj.tcga < 0.05)

