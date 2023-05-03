docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_rsession:1.3  /bin/bash # nolint: error.

#'#################################################################################
#'#################################################################################
#' Analyze PRAD
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(GSVA)
library(e1071)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(rjson)
library(NetActivity)


load("results/GSE169038/allGenes.se.RData")
# se.tcga_genes <- loadHDF5SummarizedExperiment("results/GSE169038/", prefix = "network_genes")


# genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
# path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
# path_N <- group_by(path.map, PathwayID) %>% summarize(N = n()) %>% mutate(category = PathwayID)

# prad.feat <- read.table("results/GSE169038/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
# paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
# paths.vec <- as.character(paths[, 1])
# colnames(prad.feat) <- paths.vec

# load("results/GTEx_coding/paths_filt2_full_v3.11/gtex_tcga_comparative.Rdata")
# pc.feat <- prcomp(prad.feat)

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
# prad.feat.filt <- prad.feat[!(se$primary == 1 |  se$secondary == 1), ]
# se.tcga_genes_filt <- se.tcga_genes[, !(se$primary == 1 |  se$secondary == 1)]
save(se.filt, file = "results/GSE169038/SE_filt.Rdata")

## DE genes
mod <- model.matrix(~  gleason + race + decipher, colData(se.filt))
lm.genes <- lmFit(assay(se.filt), mod) %>% eBayes()
tab.genes_geo <- topTable(lm.genes, coef = 2, n = Inf)

tab.genes_geo$gene <- rowData(se)[as.character(rownames(tab.genes_geo )), "gene"]

save(tab.genes_geo, file = "results/GSE169038/de_genes_results.Rdata")


## NetActivity
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
gse169038_prep <- prepareSummarizedExperiment(gse169038_se_filt, "gtex_gokegg")
gse169038_scores <- computeGeneSetScores(gse169038_prep, "gtex_gokegg")
save(gse169038_scores, file = "results/GSE169038/NetActivity_scores.Rdata")

mod_gse169038 <- model.matrix(~ gleason + race + decipher, colData(gse169038_scores))
lm.paths <- lmFit(assay(gse169038_scores), mod_gse169038) %>% eBayes()
tab.paths_geo <- topTable(lm.paths, coef = 2, n = Inf)
tab.paths_geo$category <- rownames(tab.paths_geo)

# tab.paths$pathway <- rownames(tab.paths)

# comb_fgsea <- left_join(fgseaRes, tab.paths, by = "pathway")
#
# png("figures/GSE169038_pval_gseavsPaths.png")
#
# ggplot(comb_fgsea, aes(x = -log10(P.Value ), y = -log10(pval ))) +
#  geom_point() +
#  scale_y_continuous(name = "-log10 P-value GSEA") +
#  scale_x_continuous(name = "-log10 P-value Pathways") +
#  theme_bw()
# dev.off()

tab.paths_geo$DE_prop <- sapply( tab.paths_geo$category, function(cat) {
  genes <- subset(path.map, PathwayID  == cat )$Symbol
  mini_tab <- subset(tab.genes_geo, gene %in% genes)
  mean(mini_tab$adj.P.Val < 0.05, na.rm = TRUE)
})

png("figures/GSE169038_propDE_vs_pvalPaths.png")
ggplot(tab.paths_geo, aes(x = DE_prop, y = -log10(P.Value ))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "-log10 P-value Pathways") +
 theme_bw()
dev.off()

png("figures/GSE169038_propDE_vs_logFCPaths.png")
ggplot(tab.paths_geo, aes(x = DE_prop , y = abs(logFC))) +
 geom_point() +
 scale_x_continuous(name = "Proportion of genes DE") +
 scale_y_continuous(name = "logFC Pathways (absolute value)") +
 theme_bw()
dev.off()


cor(tab.paths_geo$DE_prop, -log10(tab.paths_geo$P.Value))
# [1] 0.3426017

cor(abs(tab.paths_geo$logFC), tab.paths_geo$DE_prop)
# [1] 0.3162817

cor(tab.paths_geo$DE_prop[tab.paths_geo$DE_prop > 0 ], abs(tab.paths_geo$logFC)[tab.paths_geo$DE_prop > 0 ], use = "complete")
# 0.3053439



save(tab.paths_geo, file = "results/GSE169038/pathways_results.Rdata")

## GSVA
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
geo_gsva <- gsva(data.matrix(assay(se.tcga_genes_filt)), path_genes, min.sz=5, max.sz=500)

lm.gsva <- lmFit(geo_gsva, mod) %>% eBayes()
tab.gsva_geo <- topTable(lm.gsva, coef = 2, n = Inf)
tab.gsva_geo$category <- rownames(tab.gsva_geo)

save(tab.gsva_geo, geo_gsva, file = "results/GSE169038/GSVA_results.Rdata")

tab.gsva_geo_race <- topTable(lm.gsva, coef = 3, n = Inf)

tab.paths_geo_race <- topTable(lm.paths, coef = 3, n = Inf)


## hipathia
rownames(se) <- rowData(se)$gene
trans_data <- translate_data(se, "hsa")
exp_data <- normalize_data(trans_data)
hip_pathways <- load_pathways(species = "hsa")

hip.res_geo <- hipathia(exp_data, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.geo_path <- get_paths_data(hip.res_geo )
hip.comp_geo <- do_wilcoxon(hip.geo_path, hip.geo_path$gleason, g1 = "High", g2 = "Low")

save(hip.comp_geo, hip.geo_path, file = "results/GSE169038/hipathia.res.Rdata")





#
# ## Compare values between TCGA and GEO
# tcga.feat.all <- read.table("results/TCGA_gexp_coding_PRAD/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
# colnames(tcga.feat.all) <- paths.vec
#
# comb.feat.all <- rbind(tcga.feat.all, prad.feat)
# pc.comb.all <- prcomp(comb.feat.all)
#
# ## Subset data
# load("data/tcga_gexp_combat.Rdata")
# tcga.prad <- gexp_tcga_combat[as.character(genes$V1), gexp_tcga_combat$project_id == "TCGA-PRAD"]
#
#
# comb.all.df.pc <- data.frame(pc.comb.all$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat.all), nrow(prad.feat))),
#   Type = c(tcga.prad$sample_type, rep("Primary Tumor", nrow(prad.feat) ))) %>%
#   mutate(category = paste(dataset, Type))
#
# png("figures/TCGA_GEO_all_samples_PCA.png")
# ggplot(comb.all.df.pc, aes(x = PC1, y = PC2, color = category)) +
#   geom_point() +
#   theme_bw()
# dev.off()
#
# pc.tcga <- prcomp(tcga.feat)
# tcga.df.pc <- data.frame(pc.tcga$x,  Type = tcga.prad$sample_type)
# ggplot(tcga.df.pc, aes(x = PC1, y = PC2, color = Type)) +
#   geom_point() +
#   theme_bw()
#
#
# tcga.feat <- read.table("results/TCGA_gexp_coding_PRAD_tumor/comb_paths3_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
# colnames(tcga.feat) <- paths.vec
#
# tcga.prad <- tcga.prad[, !is.na(tcga.prad$paper_Reviewed_Gleason_category)]
#
# ## define TCGA gleason
# tcga.prad$gleason <- ifelse(tcga.prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")
#
# comb.feat <- rbind(tcga.feat, prad.feat.filt)
# comb.feat$dataset <- rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt)))
# comb.feat$gleason <- c(tcga.prad$gleason, se.filt$gleason)
#
#
# png("figures/TCGAvsGEO_GO0000212_boxplot.png")
# ggplot(comb.feat, aes(x = gleason, y = `GO:0000212`, color = dataset)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
#

#
#
# pc.comb <- prcomp(data.matrix(comb.feat[, 1:1337]))
#
# comb.df.pc <- data.frame(pc.comb$x, dataset = rep(c("TCGA", "GEO"), c(nrow(tcga.feat), nrow(prad.feat.filt))),
#       gleason = c(tcga.prad$gleason, se.filt$gleason))
#
# png("figures/TCGAGEO_comb_PC_dataset.png")
# ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = dataset)) +
#   geom_point() +
#   theme_bw()
# dev.off()
# #
# png("figures/TCGAGEO_comb_PC_gleason.png")
# ggplot(comb.df.pc, aes(x = PC1, y = PC2, color = gleason)) +
#   geom_point() +
#   theme_bw()
# dev.off()
#
#

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
sel_paths <- subset(tab.paths_geo , adj.P.Val < 0.05)
svm_gleason_geo_filt <- svm(gleason ~ ., df_svm[, c("gleason", gsub(":", "_", rownames( sel_paths)))])
pred.geo2_filt <- predict(svm_gleason_geo_filt, prad.feat.filt)
table(prediction = pred.geo2_filt , real = se.filt$gleason )


## GSVA
rownames(geo_gsva) <- gsub(":", "_", rownames(geo_gsva))

df_svm_gsva <-  data.frame(gleason = factor(se.filt$gleason), t(geo_gsva))
svm_gleason_geo_gsva <- svm(gleason ~ ., df_svm_gsva)
pred.geo_gsva <- predict(svm_gleason_geo_gsva, t(geo_gsva))
table(prediction = pred.geo_gsva , real = se.filt$gleason )

sel_paths2 <- subset(tab.gsva_geo  , adj.P.Val < 0.05)
svm_gleason_geo_gsva_filt <- svm(gleason ~ ., df_svm_gsva[, c("gleason", gsub(":", "_", rownames( sel_paths2)))])
pred.geo2_gsva_filt <- predict(svm_gleason_geo_gsva_filt,t(geo_gsva))
table(prediction = pred.geo2_gsva_filt , real = se.filt$gleason )

## Train SVM to classify gleason in this test

runSVM_Predicition <- function(df, train_ids, feats, true_vals){
  svm_gleason <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df[train_ids, ])
  pred <- predict(svm_gleason, feats[-train_ids, ])
  table(prediction = pred, real = true_vals[-train_ids] )
}

run_svm_iteration <- function(df_path, df_gsva, feats_path, feats_gsva, true_vals){

  train_high <- sample(which(true_vals == "High"), size = sum(true_vals == "High")*0.7)
  train_low <- sample(which(true_vals == "Low"), size = sum(true_vals == "Low")*0.7)
  train_ids <- c(train_high, train_low)

  tab_path <- runSVM_Predicition(df_path, train_ids, feats_path, true_vals)
  tab_gsva <- runSVM_Predicition(df_gsva, train_ids, feats_gsva, true_vals)
  list(path = tab_path, gsva = tab_gsva)
}

computeAccuracies <- function(tab){
  acc <- (tab[1] + tab[4])/sum(tab)
  SN <- tab[1, 1] / sum(tab[, 1])
  SP <- tab[2, 2] / sum(tab[, 2])
  c(Accuracy = acc, SN = SN, SP = SP)
}

rownames(geo_gsva) <- gsub(":", "_", rownames(geo_gsva))
colnames(prad.feat.filt) <- gsub(":", "_", colnames( prad.feat.filt))

### All features

df_svm_paths <- data.frame(gleason = factor(se.filt$gleason), prad.feat.filt)
df_svm_gsva <-  data.frame(gleason = factor(se.filt$gleason), t(geo_gsva))


set.seed(1)
iter_all <- parallel::mclapply(seq_len(50), function(i)  run_svm_iteration(df_svm_paths, df_svm_gsva, prad.feat.filt, t(geo_gsva), se.filt$gleason ), mc.cores = 10)
acc_path_all <- sapply(iter_all, function(x) computeAccuracies(x$path))
acc_gsva_all <- sapply(iter_all, function(x) computeAccuracies(x$gsva))

rowMeans(acc_path_all)
rowMeans(acc_gsva_all)

### Significant features in TCGA
tcga_paths <- gsub(":", "_", subset(comb_paths, Signif %in% c("TCGA", "Both"))$category)
tcga_gsva <- gsub(":", "_", subset(comb.gsva, Signif %in% c("TCGA", "Both"))$category)

df_svm_paths_tcga <- df_svm_paths[, c("gleason", tcga_paths)]
df_svm_gsva_tcga <- df_svm_gsva[, c("gleason", tcga_gsva)]

set.seed(1)
iter_tcga <- parallel::mclapply(seq_len(50), function(i)  run_svm_iteration(df_svm_paths_tcga, df_svm_gsva_tcga, prad.feat.filt[, tcga_paths], t(geo_gsva)[, tcga_gsva], se.filt$gleason ), mc.cores = 10)
acc_path_tcga <- sapply(iter_tcga, function(x) computeAccuracies(x$path))
acc_gsva_tcga <- sapply(iter_tcga, function(x) computeAccuracies(x$gsva))

rowMeans(acc_path_tcga)
rowMeans(acc_gsva_tcga)

class_df <- rbind(data.frame(t(acc_path_all), Features = "All", Model = "GO + KEGG model"),
                  data.frame(t(acc_gsva_all), Features = "All", Model = "GSVA"),
                  data.frame(t(acc_path_tcga), Features = "TCGA", Model = "GO + KEGG model"),
                  data.frame(t(acc_gsva_tcga), Features = "TCGA", Model = "GSVA")) %>%
                  as_tibble()



png("figures/GEO_classification.png")
class_df %>%
  gather(Measure, Value, 1:3) %>%
  ggplot(class_df, aes(x = Measure, y = Value, color = Model)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(~ Features)
dev.off()

summary(lm(SN ~ Model, class_df, subset = Features == "All"))
summary(lm(SP ~ Model, class_df, subset = Features == "All"))
summary(lm(Accuracy ~ Model, class_df, subset = Features == "All"))

summary(lm(SN ~ Model, class_df, subset = Features == "TCGA"))
summary(lm(SP ~ Model, class_df, subset = Features == "TCGA"))
summary(lm(Accuracy ~ Model, class_df, subset = Features == "TCGA"))

### Gene sets with similar training in both
sel_paths <- gsub(":", "_", subset(path_df, abs(cor_path) > 0.7)$path)

df_svm_paths_sel <- df_svm_paths[, c("gleason", sel_paths)]
df_svm_gsva_sel <- df_svm_gsva[, c("gleason", sel_paths)]

set.seed(1)
iter_sel <- parallel::mclapply(seq_len(50), function(i)  run_svm_iteration(df_svm_paths_sel, df_svm_gsva_sel, prad.feat.filt[, sel_paths], t(geo_gsva)[, sel_paths], se.filt$gleason ), mc.cores = 10)
acc_path_sel <- sapply(iter_sel, function(x) computeAccuracies(x$path))
acc_gsva_sel <- sapply(iter_sel, function(x) computeAccuracies(x$gsva))

rowMeans(acc_path_sel)
rowMeans(acc_gsva_sel)





train_high <- sample(which(se.filt$gleason == "High"), size = sum(se.filt$gleason == "High")*0.7)
train_low <- sample(which(se.filt$gleason == "Low"), size = sum(se.filt$gleason == "Low")*0.7)

df_svm <- data.frame(gleason = factor(se.filt$gleason), prad.feat.filt)
svm_gleason_geo <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df_svm[c(train_high, train_low), ])
pred.geo2 <- predict(svm_gleason_geo, prad.feat.filt[-c(train_high, train_low), ])
table(prediction = pred.geo2 , real = se.filt[, -c(train_high, train_low)]$gleason )

df_svm_gsva <-  data.frame(gleason = factor(se.filt$gleason), t(geo_gsva))
svm_gleason_geo_gsva <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df_svm_gsva[c(train_high, train_low), ])
pred.geo_gsva <- predict(svm_gleason_geo_gsva, t(geo_gsva)[-c(train_high, train_low), ])
table(prediction = pred.geo_gsva , real = se.filt[, -c(train_high, train_low)]$gleason )

## Selected features
sel_paths <- subset(tab.paths_tcga , adj.P.Val < 0.05)
svm_gleason_geo_filt <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df_svm[c(train_high, train_low), c("gleason", gsub(":", "_", rownames( sel_paths)))])
pred.geo2_filt <- predict(svm_gleason_geo_filt, prad.feat.filt[-c(train_high, train_low), ])
table(prediction = pred.geo2_filt , real = se.filt[, -c(train_high, train_low)]$gleason )

sel_paths2 <- subset(tab.gsva_tcga, adj.P.Val < 0.05)
svm_gleason_geo_gsva_filt <- svm(gleason ~ .,  class.weights = c(Low = 0.1, High = 0.9), df_svm_gsva[c(train_high, train_low), c("gleason", gsub(":", "_", rownames( sel_paths2)))])
pred.geo2_gsva_filt <- predict(svm_gleason_geo_gsva_filt,t(geo_gsva)[-c(train_high, train_low), ])
table(prediction = pred.geo2_gsva_filt , real = se.filt[, -c(train_high, train_low)]$gleason  )

save(svm_gleason_geo, svm_gleason_geo_filt, svm_gleason_geo_gsva_filt, svm_gleason_geo_gsva, file = "results/TCGA_PRAD/svm_model_geo.Rdata")
