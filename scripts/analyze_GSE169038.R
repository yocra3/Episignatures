#'#################################################################################
#'#################################################################################
#' Analyze GSE169038 datasete
#'#################################################################################
#'#################################################################################


## Load libraries
library(fgsea)
library(limma)
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)
library(e1071)

load("results/GSE169038/allGenes.se.RData")
se.tcga_genes <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")

load( "results/TCGA_PRAD/pathways_results.Rdata")

genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
path.map <- read.table("results/preprocess/go_kegg_final_gene_map.tsv", header = TRUE)

prad.feat <- read.table("results/GSE169038/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(prad.feat) <- paths.vec


pc.feat <- prcomp(prad.feat)

se$race <- ifelse(grepl( "White", se$characteristics_ch1.4), "EUR", "AFR")
se$decipher <- factor(gsub("decipher risk group: ", "", se$characteristics_ch1.3), levels = c("Lower", "Average", "Higher"))
se$primary <- gsub("primary gleason: ", "", se$characteristics_ch1.1)
se$secondary <- gsub("secondary gleason: ", "", se$characteristics_ch1.2)
se$primary <- as.numeric(ifelse(se$primary == "--", 1, se$primary))
se$secondary <- as.numeric(ifelse(se$secondary == "--", 1, se$secondary))
se$gleason_cat <- paste(se$primary, se$secondary, sep = "-")
se$gleason <- ifelse(se$primary == 5 |  se$secondary == 5 | se$gleason_cat == "4+4", "High", "Low")

## Subset samples with gleason < 3
se.filt <- se[, !(se$primary == 1 |  se$secondary == 1)]
prad.feat.filt <- prad.feat[!(se$primary == 1 |  se$secondary == 1), ]


## DE genes
mod <- model.matrix(~  gleason + race + decipher, colData(se.filt))
lm.genes <- lmFit(assay(se.filt), mod) %>% eBayes()
tab.genes <- topTable(lm.genes, coef = 2, n = Inf)


## fgsea
ranks <- tab.genes$logFC
names(ranks) <- rowData(se)[rownames(tab.genes), "gene"]

names(paths.vec) <- paths.vec
pathways <- lapply( paths.vec, function(x) subset(path.map, PathwayID  == x )$Symbol)
fgseaRes <- fgsea(pathways = pathways, stats = ranks)
# fgseaRes$category <- fgseaRes$pathway


## DE paths
lm.paths <- lmFit(t(prad.feat.filt), mod) %>% eBayes()
tab.paths <- topTable(lm.paths, coef = 2, n = Inf)
tab.paths$category <- rownames(tab.paths)
tab.paths$pathway <- rownames(tab.paths)

comb_fgsea <- left_join(fgseaRes, tab.paths, by = "pathway")

png("figures/GSE169038_pval_gseavsPaths.png")

ggplot(comb_fgsea, aes(x = -log10(P.Value ), y = -log10(pval ))) +
 geom_point() +
 scale_y_continuous(name = "-log10 P-value GSEA") +
 scale_x_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()

tab.genes$gene <- rowData(se)[rownames(tab.genes), "gene"]

tab.paths$DE_prop <- sapply( tab.paths$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- subset(tab.genes, gene %in% genes)
  mean(mini_tab$adj.P.Val < 0.05)
})

png("figures/GSE169038_propDE_vs_pvalPaths.png")
ggplot(tab.paths, aes(x = DE_prop, y = -log10(P.Value ))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()


cor(tab.paths$ DE_prop, -log10(tab.paths$P.Value))
# 0.2668079

load("results/TCGA_PRAD/pathways_results.Rdata")
comb_paths <- left_join(tab.path, tab.paths, by = "category", suffix = c(".TCGA", ".GEO")) %>%
  mutate(Signif = ifelse(adj.P.Val.TCGA < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))


png("figures/TCGAvsGEO_logFC.png")
ggplot(comb_paths, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw()
dev.off()

png("figures/TCGAvsGEO_pval.png")
ggplot(comb_paths, aes(x = -log10(P.Value.TCGA), y =  -log10(P.Value.GEO), col = Signif)) +
  geom_point() +
  theme_bw()
dev.off()

summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths ))
summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths, subset = Signif != "None" ))

cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)
cor(comb_paths[comb_paths$Signif != "None", ]$logFC.TCGA, comb_paths[comb_paths$Signif != "None", ]$logFC.GEO)
table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif), margin = 2)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif != "None"), margin = 2)

## Compare values between TCGA and GEO
tcga.feat.all <- read.table("results/TCGA_gexp_coding_PRAD/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(tcga.feat.all) <- paths.vec

comb.feat.all <- rbind(tcga.feat.all, prad.feat)
pc.comb.all <- prcomp(comb.feat.all)

## Subset data
load("data/tcga_gexp_combat.Rdata")
tcga.prad <- gexp_tcga_combat[as.character(genes$V1), gexp_tcga_combat$project_id == "TCGA-PRAD"]


comb.all.df.pc <- data.frame(pc.comb.all$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat.all), nrow(prad.feat))),
  Type = c(tcga.prad$sample_type, rep("Primary Tumor", nrow(prad.feat) ))) %>%
  mutate(category = paste(dataset, Type))

png("figures/TCGA_GEO_all_samples_PCA.png")
ggplot(comb.all.df.pc, aes(x = PC1, y = PC2, color = category)) +
  geom_point() +
  theme_bw()
dev.off()

pc.tcga <- prcomp(tcga.feat)
tcga.df.pc <- data.frame(pc.tcga$x,  Type = tcga.prad$sample_type)
ggplot(tcga.df.pc, aes(x = PC1, y = PC2, color = Type)) +
  geom_point() +
  theme_bw()


tcga.feat <- read.table("results/TCGA_gexp_coding_PRAD_tumor/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(tcga.feat) <- paths.vec

tcga.prad <- tcga.prad[, !is.na(tcga.prad$paper_Reviewed_Gleason_category)]

## define TCGA gleason
tcga.prad$gleason <- ifelse(tcga.prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")

comb.feat <- rbind(tcga.feat, prad.feat.filt)
comb.feat$dataset <- rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt)))
comb.feat$gleason <- c(tcga.prad$gleason, se.filt$gleason)


png("figures/TCGAvsGEO_GO0000212_boxplot.png")
ggplot(comb.feat, aes(x = gleason, y = `GO:0000212`, color = dataset)) +
  geom_boxplot() +
  theme_bw()
dev.off()




pc.comb <- prcomp(data.matrix(comb.feat[, 1:1337]))

comb.df.pc <- data.frame(pc.comb$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt))),
      gleason = c(tcga.prad$gleason, se.filt$gleason))

png("figures/TCGAGEO_comb_PC_dataset.png")
ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point() +
  theme_bw()
dev.off()
#
png("figures/TCGAGEO_comb_PC_gleason.png")
ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = gleason)) +
  geom_point() +
  theme_bw()
dev.off()



## Test SVM to classify gleason
load("results/TCGA_PRAD/svm_model.Rdata")
colnames(prad.feat.filt) <- gsub(":", "_", colnames( prad.feat.filt))
pred.geo <- predict(svm_gleason, prad.feat.filt)
table(prediction = pred.geo , real = se.filt$gleason )

pred.geo_filt <- predict(svm_gleason_filt, prad.feat.filt)
table(prediction = pred.geo_filt , real = se.filt$gleason )


## Train SVM to classify gleason
### All features
df_svm <-  data.frame(gleason = factor(se.filt$gleason), prad.feat.filt)
svm_gleason_geo <- svm(gleason ~ ., df_svm)
pred.geo2 <- predict(svm_gleason_geo, prad.feat.filt)
table(prediction = pred.geo2 , real = se.filt$gleason )

## Selected features
sel_paths <- subset(tab.paths, adj.P.Val < 0.05)
svm_gleason_geo_filt <- svm(gleason ~ ., df_svm[, c("gleason", gsub(":", "_", rownames( sel_paths)))])
pred.geo2_filt <- predict(svm_gleason_geo_filt, prad.feat.filt)
table(prediction = pred.geo2_filt , real = se.filt$gleason )

save(svm_gleason_geo, svm_gleason_geo_filt, file = "results/TCGA_PRAD/svm_model_geo.Rdata")


# Compare Gene DE
load("results/TCGA_PRAD/genes_results.Rdata")
res$gene <- rownames(res)
tab.genes$gene <- rowData(se)[rownames(tab.genes), "gene"]
comb.genes <- left_join(data.frame(res), tab.genes, by = "gene")   %>%
  mutate(Signif = ifelse(!is.na(padj) & padj  < 0.05, ifelse(adj.P.Val < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val < 0.05, "GEO", "None"))) %>%
  filter(!is.na(pvalue ) & !is.na(P.Value  ))

png("figures/TCGAvsGEO_genes_logFC.png")
ggplot(comb.genes, aes(x = log2FoldChange, y = logFC, col = Signif)) +
  geom_point() +
  theme_bw()
dev.off()


cor(comb.genes$log2FoldChange, comb.genes$logFC)
# [1] 0.2488376
cor(comb.genes[comb.genes$Signif != "None", ]$log2FoldChange, comb.genes[comb.genes$Signif != "None", ]$logFC, use = "complete")
# [1] 0.3250043

table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif), margin = 2)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif != "None"), margin = 2)
