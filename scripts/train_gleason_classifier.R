#'#################################################################################
#'#################################################################################
#' Train gleason
#'#################################################################################
#'#################################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(e1071)
library(HDF5Array)
library(org.Hs.eg.db)
library(parallel)
library(rjson)


## Load GSE169038
load("results/GSE169038/allGenes.se.RData")

genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)

geo.feat <- read.table("results/GSE169038/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(geo.feat) <- paths.vec

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
geo.feat.filt <- geo.feat[!(se$primary == 1 |  se$secondary == 1), ]
hip.geo_path.filt <- hip.geo_path[, colnames(se.filt)]
load("results/GSE169038/GSVA_results.Rdata")
load( "results/GSE169038/pathways_results.Rdata")

## Load PRAD
load("data/tcga_gexp_combat.Rdata")
prad.feat <- read.table("results/GTEx_coding_PRAD_tumor/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(prad.feat) <- paths.vec

prad <- gexp_tcga_combat[genes$V1, gexp_tcga_combat$project_id == "TCGA-PRAD"]
prad <- prad[, !is.na(prad$paper_Reviewed_Gleason_category)]
prad$gleason <- ifelse(prad$paper_Reviewed_Gleason_category == ">=8", "High", "Low")

load("results/TCGA_PRAD/pathways_results.Rdata")
load("results/TCGA_PRAD/GSVA_results.Rdata")
load("results/GSE169038/hipathia.res.Rdata")
load("results/TCGA_PRAD/hipathia.res.Rdata")


## Combine results
comb_paths <- left_join(tab.path_prad, tab.paths_geo, by = "category", suffix = c(".TCGA", ".GEO")) %>%
  as_tibble() %>%
  mutate(Signif = ifelse(adj.P.Val.TCGA < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))


comb.gsva <- inner_join(tab.gsva_prad, tab.gsva_geo, by = "category", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(adj.P.Val.TCGA) & adj.P.Val.TCGA  < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))

#
comb.hipathia <- inner_join(hip.comp_prad, hip.comp_geo, by = "name", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(FDRp.value.TCGA) & FDRp.value.TCGA  < 0.05, ifelse(FDRp.value.GEO < 0.05, "Both", "TCGA"),
                              ifelse(FDRp.value.GEO < 0.05, "GEO", "None"))) %>%
   left_join(mutate(data.frame(rowData(hip.geo_path)), name = feat.name), by = "name")



## Train SVM to classify gleason in this test
runSVM_Predicition <- function(df, train_ids, feats, true_vals, prop = c(Low = 0.1, High = 0.9)){
  svm_gleason <- svm(gleason ~ ., class.weights = prop, df[train_ids, ])
  pred <- predict(svm_gleason, feats[-train_ids, ])
  table(prediction = pred, real = true_vals[-train_ids] )
}

run_svm_iteration <- function(df_path, df_gsva, df_hip, feats_path, feats_gsva, feats_hip, true_vals, prop = c(Low = 0.1, High = 0.9)){

  train_high <- sample(which(true_vals == "High"), size = sum(true_vals == "High")*0.7)
  train_low <- sample(which(true_vals == "Low"), size = sum(true_vals == "Low")*0.7)
  train_ids <- c(train_high, train_low)

  tab_path <- runSVM_Predicition(df_path, train_ids, feats_path, true_vals, prop)
  tab_gsva <- runSVM_Predicition(df_gsva, train_ids, feats_gsva, true_vals, prop)
  tab_hip <- runSVM_Predicition(df_hip, train_ids, feats_hip, true_vals, prop)

  list(path = tab_path, gsva = tab_gsva, hip = tab_hip)
}

computeAccuracies <- function(tab){
  acc <- (tab[1] + tab[4])/sum(tab)
  SN <- tab[1, 1] / sum(tab[, 1])
  SP <- tab[2, 2] / sum(tab[, 2])
  c(Accuracy = acc, SN = SN, SP = SP)
}

rownames(geo_gsva) <- gsub(":", "_", rownames(geo_gsva))
rownames(prad_gsva) <- gsub(":", "_", rownames(prad_gsva))

colnames(geo.feat.filt) <- gsub(":", "_", colnames( geo.feat.filt))
colnames(prad.feat) <- gsub(":", "_", colnames( prad.feat))

hip.geo_path.filt <- hip.geo_path[, colnames(geo_gsva)]

rownames(hip.geo_path.filt) <- gsub("-", "_", rownames(hip.geo_path.filt))
rownames(hip.geo_path.filt) <- gsub(" ", "_", rownames(hip.geo_path.filt))

rownames(hip.prad_vals) <- gsub("-", "_", rownames(hip.prad_vals ))
rownames(hip.prad_vals) <- gsub(" ", "_", rownames(hip.prad_vals))

## Train in GEO
### All features

df_svm_paths <- data.frame(gleason = factor(se.filt$gleason), geo.feat.filt)
df_svm_gsva <-  data.frame(gleason = factor(se.filt$gleason), t(geo_gsva))
df_svm_hip <-  data.frame(gleason = factor(se.filt$gleason), t(assay(hip.geo_path.filt)))


set.seed(1)
iter_all <- parallel::mclapply(seq_len(50), function(i)
  run_svm_iteration(df_svm_paths, df_svm_gsva, df_svm_hip, geo.feat.filt, t(geo_gsva), t(assay(hip.geo_path.filt)), se.filt$gleason ), mc.cores = 10)
acc_path_all <- sapply(iter_all, function(x) computeAccuracies(x$path))
acc_gsva_all <- sapply(iter_all, function(x) computeAccuracies(x$gsva))
acc_hip_all <- sapply(iter_all, function(x) computeAccuracies(x$hip))

rowMeans(acc_path_all)
rowMeans(acc_gsva_all)
rowMeans(acc_hip_all)

### Significant features in TCGA
tcga_paths <- gsub(":", "_", subset(comb_paths, Signif %in% c("TCGA", "Both"))$category)
tcga_gsva <- gsub(":", "_", subset(comb.gsva, Signif %in% c("TCGA", "Both"))$category)
tcga_hip <- gsub("-", "_", subset(comb.hipathia, Signif %in% c("TCGA", "Both"))$feat.ID)
tcga_hip <- gsub(" ", "_", tcga_hip)

df_svm_paths_tcga <- df_svm_paths[, c("gleason", tcga_paths)]
df_svm_gsva_tcga <- df_svm_gsva[, c("gleason", tcga_gsva)]
df_svm_hip_tcga <- df_svm_hip[, c("gleason", tcga_hip)]

set.seed(1)
iter_tcga <- parallel::mclapply(seq_len(50), function(i)
  run_svm_iteration(df_svm_paths_tcga, df_svm_gsva_tcga, df_svm_hip_tcga, geo.feat.filt[, tcga_paths], t(geo_gsva)[, tcga_gsva], t(assay(hip.geo_path.filt[tcga_hip, ])), se.filt$gleason ), mc.cores = 10)
acc_path_tcga <- sapply(iter_tcga, function(x) computeAccuracies(x$path))
acc_gsva_tcga <- sapply(iter_tcga, function(x) computeAccuracies(x$gsva))
acc_hip_tcga <- sapply(iter_tcga, function(x) computeAccuracies(x$hip))

rowMeans(acc_path_tcga)
rowMeans(acc_gsva_tcga)
rowMeans(acc_hip_tcga)

class_df <- rbind(data.frame(t(acc_path_all), Features = "All", Model = "GO + KEGG model"),
                  data.frame(t(acc_gsva_all), Features = "All", Model = "GSVA"),
                  data.frame(t(acc_hip_tcga), Features = "All", Model = "Hipathia"),
                  data.frame(t(acc_path_tcga), Features = "Selected", Model = "GO + KEGG model"),
                  data.frame(t(acc_hip_tcga), Features = "Selected", Model = "Hipathia"),
                  data.frame(t(acc_gsva_tcga), Features = "Selected", Model = "GSVA")) %>%
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


## Train in TCGA
### All features
df_svm_paths2 <- data.frame(gleason = factor(prad$gleason), prad.feat)
df_svm_gsva2 <-  data.frame(gleason = factor(prad$gleason), t(assay(prad_gsva)))
df_svm_hip2 <-  data.frame(gleason = factor(prad$gleason), t(assay(hip.prad_vals)))

set.seed(1)
iter_all2 <- parallel::mclapply(seq_len(50), function(i)
  run_svm_iteration(df_svm_paths2, df_svm_gsva2, df_svm_hip2, prad.feat, t(assay(prad_gsva)), t(assay(hip.prad_vals)), prad$gleason, prop = c(High = 0.75, Low = 0.25) ), mc.cores = 10)
acc_path_all2 <- sapply(iter_all2, function(x) computeAccuracies(x$path))
acc_gsva_all2 <- sapply(iter_all2, function(x) computeAccuracies(x$gsva))
acc_hip_all2 <- sapply(iter_all2, function(x) computeAccuracies(x$hip))

rowMeans(acc_path_all2)
rowMeans(acc_gsva_all2)
rowMeans(acc_hip_all2)

### Significant features in TCGA
tcga_paths2 <- gsub(":", "_", subset(comb_paths, Signif %in% c("GEO", "Both"))$category)
tcga_gsva2 <- gsub(":", "_", subset(comb.gsva, Signif %in% c("GEO", "Both"))$category)
tcga_hip2 <- gsub("-", "_", subset(comb.hipathia, Signif %in% c("GEO", "Both"))$feat.ID)
tcga_hip2 <- gsub(" ", "_", tcga_hip2)

df_svm_paths_tcga2 <- df_svm_paths2[, c("gleason", tcga_paths2)]
df_svm_gsva_tcga2 <- df_svm_gsva2[, c("gleason", tcga_gsva2)]
df_svm_hip_tcga2 <- df_svm_hip2[, c("gleason", tcga_hip2)]

set.seed(1)
iter_tcga2 <- parallel::mclapply(seq_len(50), function(i)
  run_svm_iteration(df_svm_paths_tcga2, df_svm_gsva_tcga2, df_svm_hip_tcga2, prad.feat[, tcga_paths2], t(assay(prad_gsva))[, tcga_gsva2], t(assay(hip.prad_vals))[, tcga_hip2], prad$gleason, prop = c(High = 0.75, Low = 0.25)), mc.cores = 10)
acc_path_tcga2 <- sapply(iter_tcga2, function(x) computeAccuracies(x$path))
acc_gsva_tcga2 <- sapply(iter_tcga2, function(x) computeAccuracies(x$gsva))
acc_hip_tcga2 <- sapply(iter_tcga2, function(x) computeAccuracies(x$hip))

rowMeans(acc_path_tcga2)
rowMeans(acc_gsva_tcga2)
rowMeans(acc_hip_tcga2)

class_df2 <- rbind(data.frame(t(acc_path_all2), Features = "All", Model = "GO + KEGG model"),
                  data.frame(t(acc_gsva_all2), Features = "All", Model = "GSVA"),
                  data.frame(t(acc_hip_tcga2), Features = "All", Model = "Hipathia"),
                  data.frame(t(acc_path_tcga2), Features = "Selected", Model = "GO + KEGG model"),
                  data.frame(t(acc_gsva_tcga2), Features = "Selected", Model = "GSVA"),
                  data.frame(t(acc_hip_tcga2), Features = "Selected", Model = "Hipathia")) %>%
                  as_tibble()



png("figures/TCGA_classification.png")
class_df2 %>%
  gather(Measure, Value, 1:3) %>%
  ggplot(aes(x = Measure, y = Value, color = Model)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(~ Features)
dev.off()


png("figures/gleason_classification.png", width = 600)
rbind(mutate(class_df, Dataset = "GEO"),
      mutate(class_df2, Dataset = "TCGA")) %>%
  mutate(Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
  gather(Measure, Value, 1:3) %>%
  ggplot(aes(x = Measure, y = Value, color = Model)) +
    geom_boxplot() +
    theme_bw() +
    facet_grid(Features ~ Dataset)
dev.off()



summary(lm(SN ~ Model, class_df2, subset = Features == "All"))
summary(lm(SP ~ Model, class_df2, subset = Features == "All"))
summary(lm(Accuracy ~ Model, class_df2, subset = Features == "All"))

summary(lm(SN ~ Model, class_df2, subset = Features == "GEO"))
summary(lm(SP ~ Model, class_df2, subset = Features == "GEO"))
summary(lm(Accuracy ~ Model, class_df2, subset = Features == "GEO"))

summary(lm(Accuracy ~ Features, class_df, subset = Model == "GO + KEGG model"))
summary(lm(Accuracy ~ Features, class_df2, subset = Model == "GO + KEGG model"))

summary(lm(SN ~ Features, class_df, subset = Model == "GO + KEGG model"))
summary(lm(SN ~ Features, class_df2, subset = Model == "GO + KEGG model"))

summary(lm(SP ~ Features, class_df, subset = Model == "GO + KEGG model"))
summary(lm(SP ~ Features, class_df2, subset = Model == "GO + KEGG model"))



### Check training in one dataset and predict in other
svm_geo_paths <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df_svm_paths)
table(prediction = predict(svm_geo_paths, prad.feat), real = prad$gleason )

svm_geo_gsva <- svm(gleason ~ ., class.weights = c(Low = 0.1, High = 0.9), df_svm_gsva)
table(prediction = predict(svm_geo_gsva, t(assay(prad_gsva))), real = prad$gleason )

svm_tcga_paths <- svm(gleason ~ ., class.weights = c(Low = 0.25, High = 0.75), df_svm_paths2)
table(prediction = predict(svm_tcga_paths, geo.feat.filt), real = se.filt$gleason )

svm_tcga_gsva <- svm(gleason ~ ., class.weights = c(Low = 0.25, High = 0.75), df_svm_gsva2)
table(prediction = predict(svm_tcga_gsva, t(geo_gsva)), real =  se.filt$gleason )
