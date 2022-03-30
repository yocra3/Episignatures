#'#################################################################################
#'#################################################################################
#' Explore GSE57945 features from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(pheatmap)
library(caTools)
library(e1071)
library(hipathia)
library(rjson)


makePCdf <- function(seobj, vars){
  pc <- prcomp(t(assay(seobj)))
  pcdf <- data.frame(pc$x[, 1:10])
  pcdf <- cbind(pcdf, colData(seobj)[, vars])
}
makePCplot <- function(pcdf, var){
  ggplot(pcdf, aes_string(x = "PC1", y = "PC2", col = var)) +
  geom_point() +
  theme_bw()
}
readFeatures <- function(path, seobj){
  tab <- read.table(path, header = TRUE)
  se <- SummarizedExperiment(assay = t(data.matrix(tab)),
                              colData = colData(seobj))
}

makeHeatmap <- function(seObj){
  col_colors <- list(
    sample_type = c("Primary Tumor" = "lightgreen", "Solid Tissue Normal" = "black")
    # sex = c("Female" = "purple", "Male" = "lightblue")
  )

  pheatmap(assay(seObj), scale = "row",
           annotation_col  = data.frame(colData(seObj)[, c("sample_type"), drop = FALSE]),
           annotation_colors =  col_colors,
          show_colnames = FALSE)

}
trainSVM <-  function(seobj){
  mat <- data.matrix(assay(seobj))
  df <- data.frame(pathClass = factor(seobj$diagnosis2), t(mat))
  model_svm <- svm(pathClass ~ ., df)
}


computeTPR <- function(tumor, se){

  sub.se <- se[, se$project_id == tumor & se$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  mod <- model.matrix(~sample_type, colData(sub.se))
  paths.lm <- lmFit(assay(sub.se), mod)
  paths.lm <- eBayes(paths.lm)
  paths.top <- topTable(paths.lm, n = Inf, coef = 2)
  sum(paths.top[cancer.kegg, ]$adj.P.Val < 0.05)


}


computeLimma <- function(tumor, se){

  sub.se <- se[, se$project_id == tumor & se$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  mod <- model.matrix(~sample_type, colData(sub.se))
  paths.lm <- lmFit(t(scale(t(assay(sub.se)))), mod)
  paths.lm <- eBayes(paths.lm)
  paths.top <- topTable(paths.lm, n = Inf, coef = 2)
  paths.top$pathID <- rownames(paths.top)
  paths.top$project <- tumor
  paths.top

}


computeLimmaFPR <- function(tumor, se){

  sub.se <- se[, se$project_id == tumor & se$sample_type %in% c("Solid Tissue Normal")]
  sub.se$var <- sample(c("A", "B"), ncol(sub.se), replace = TRUE)
  mod <- model.matrix(~var, colData(sub.se))
  paths.lm <- lmFit(t(scale(t(assay(sub.se)))), mod)
  paths.lm <- eBayes(paths.lm)
  paths.top <- topTable(paths.lm, n = Inf, coef = 2)
  paths.top$pathID <- rownames(paths.top)
  paths.top$project <- tumor
  paths.top

}


getTumorPaths <- function(tumor, se){

  sub.se <- se[, se$project_id == tumor & se$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  mod <- model.matrix(~sample_type, colData(sub.se))
  paths.lm <- lmFit(assay(sub.se), mod)
  paths.lm <- eBayes(paths.lm)
  paths.top <- topTable(paths.lm, n = Inf, coef = 2)
  paths.top[cancer.kegg, ]
}

computeTPRhipathia <- function(tumor, se){

  sub.se <- se[, se$project_id == tumor & se$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]

  sample_group <- sub.se$sample_type
  comp <- do_wilcoxon(sub.se, sample_group, g1 = "Primary Tumor", g2 = "Solid Tissue Normal")
  pathways_summary <- get_pathways_summary(comp, pathways)
  subset(pathways_summary, id_pathways %in% cancer.keg)$percent_significant_paths

}

read_training <- function(path, name){

  mat <- read.table(path, row.names = 1) %>% t()
  df <- mat[, c(1, 3)] %>%
      data.frame() %>%
      mutate(epoch = seq_len(nrow(mat))) %>%
      gather(measure, mse, 1:2) %>%
      mutate(model = name)
}

## Compare trainings
## No primed values
# raw.vals <- read_training("results/TCGA_gexp_norm/kegg_v1.3/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - base")
# drop20.vals <- read_training("results/TCGA_gexp_norm/kegg_v2.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - dropout 20%")
# drop50.vals <- read_training("results/TCGA_gexp_norm/kegg_v2.2/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - dropout 50%")
# all.vals <- Reduce(rbind, list(auto.vals, raw.vals, drop20.vals, drop50.vals, primed.vals,hipathia.vals))

auto.base <- read_training("results/TCGA_gexp_combat_coding_std/autoencod_v2.2/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Pathway")
auto.post <- read_training("results/TCGA_gexp_combat_coding_std/autoencod_v3.2/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Pathway + Dense")
auto.pre <- read_training("results/TCGA_gexp_combat_coding_std/autoencod_v3.1/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Dense + Pathway")
auto.pre_post <- read_training("results/TCGA_gexp_combat_coding_std/autoencod_v4.1/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Dense + Pathway + Dense")

base <- read_training("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Pathway")
post <- read_training("results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Pathway + Dense")
pre <- read_training("results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Dense + Pathway")
pre_post <- read_training("results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_trained/TCGA_gexp_combat_coding_std_training_evaluation.txt", "Dense + Pathway + Dense")

all.vals <- rbind(Reduce(rbind, list(base, post, pre, pre_post)) %>%
  mutate(model = factor(model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense" ))) %>%
  mutate(type = "KEGG pathways"),
  Reduce(rbind, list(auto.base, auto.post, auto.pre, auto.pre_post)) %>%
    mutate(model = factor(model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense" ))) %>%
    mutate(type = "Autoencoder"))



png("figures/TCGA.kegg.training.png", width = 900)
ggplot(all.vals, aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_grid(type ~ model) +
  theme_bw()
dev.off()


png("figures/TCGA.kegg.trainingb.png", width = 600)
ggplot(all.vals, aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_grid(model ~ type) +
  theme_bw()
dev.off()


png("figures/TCGA.kegg.training.epoch10.png", width = 900)
filter(all.vals, epoch > 10) %>%
ggplot(aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_grid(type ~ model) +
  theme_bw()
dev.off()

png("figures/TCGA.kegg.training.epoch10b.png", width = 600)
filter(all.vals, epoch > 10) %>%
ggplot(aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_grid(model ~ type) +
  theme_bw()
dev.off()

## Hipathia structure
# hipathia.vals <- read_training("results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "hipathia - primed")


kegg.annot <- fromJSON(file = "data/kegg_pathways.json")
kegg.df <- lapply(kegg.annot$children, function(x) {
  top_cat <- x$name
  paths.df <- lapply(x$children, function(y){
    cat2 <- y$name
    paths <- sapply(y$children, function(z) z$name)
    data.frame(top_cat = top_cat, category = cat2, path = paths)
  })
  Reduce(rbind, paths.df) %>%
    as_tibble()
}) %>%
 Reduce(rbind, .)
kegg.df <- kegg.df %>%
  mutate(pathID = paste0("path:hsa", substring(path, 0, 5)),
          pathName = gsub("^[0-9]*  ", "", path))

## Load vsd data ####
vsd.ori <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

cancer.keg <- c("hsa04010", "hsa04310", "hsa04350", "hsa04370", "hsa04630", "hsa04024", "hsa04151", "hsa04150", "hsa04110", "hsa04210", "hsa04115", "hsa04510", "hsa04520", "hsa03320")
cancer.kegg <- paste0("path:", cancer.keg)
tumors <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA", "TCGA-UCEC")

paths <- read.table("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
cancer.kegg <- cancer.kegg[cancer.kegg %in% paths.vec]

## Pathways
se.base <- readFeatures("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.base) <- as.character(paths[, 1])
tpr_base <- sapply(tumors, computeTPR, se =  se.base )
limma_base <- lapply(tumors, computeLimma, se =  se.base )
limma_base_join <- Reduce(rbind, limma_base) %>%
  left_join(kegg.df, by = "pathID")

a <- select(limma_base_join, pathName, logFC, project) %>% spread(project, logFC)
mat <- a[, -1] %>% data.matrix()
rownames(mat) <- a$pathName

kegg.df2 <- kegg.df  %>% data.frame()
rownames(kegg.df2) <- kegg.df2$pathName


png("figures/TCGA.base.pathways_assoc_heatmap.png", width = 2500, height = 1500)
pheatmap(t(mat), annotation_col = kegg.df2[rownames(mat), c("top_cat", "category")], scale = "none")
dev.off()


## Pathways + Dense
se.post <- readFeatures("results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.post) <- as.character(paths[, 1])
tpr_post <- sapply(tumors, computeTPR, se =  se.post )

## Dense + Pathways
se.pre <- readFeatures("results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_features/prune_low_magnitude_dense_1.tsv", vsd.ori)
rownames(se.pre) <- as.character(paths[, 1])
tpr_pre <- sapply(tumors, computeTPR, se =  se.pre )
limma_pre <- lapply(tumors, computeLimma, se =  se.pre )
limma_pre_join <- Reduce(rbind, limma_pre) %>%
  left_join(kegg.df, by = "pathID")

a <- select(limma_pre_join, pathName, logFC, project) %>% spread(project, logFC)
mat <- a[, -1] %>% data.matrix()
rownames(mat) <- a$pathName

kegg.df2 <- kegg.df  %>% data.frame()
rownames(kegg.df2) <- kegg.df2$pathName


png("figures/TCGA.pathways_assoc_heatmap.png", width = 2500, height = 1500)
pheatmap(t(mat), annotation_col = kegg.df2[rownames(mat), c("top_cat", "category")], scale = "none")
dev.off()

limma_fpr <- lapply(tumors, computeLimmaFPR, se =  se.pre )
limma_fpr_join <- Reduce(rbind, limma_fpr) %>%
  left_join(kegg.df, by = "pathID")



## Dense + Pathways + Dense
se.pre_post <- readFeatures("results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_features/prune_low_magnitude_dense_1.tsv", vsd.ori)
rownames(se.pre_post) <- as.character(paths[, 1])
tpr_pre_post <- sapply(tumors, computeTPR, se =  se.pre_post )



## Plots
png("figures/TCGA.kegg_comb.tpr.png")
Reduce(rbind, list(tpr_base, tpr_post, tpr_pre, tpr_pre_post)) %>%
  as_tibble() %>%
  mutate(Model = factor(c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense"),
                          levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense"))) %>%
  gather(tumor, pos, seq_len(length(tumors))) %>%
  mutate(TPR = pos/length(cancer.kegg)) %>%
    ggplot(aes(x = Model, fill = factor(TPR))) +
    geom_bar() +
    theme_bw()
dev.off()



## From here, old!!
##################################################################################################
paths <- read.table("results/TCGA_gexp_kegg/v1.1/model_trained/pathways_names.txt", header = TRUE)


## Load autoencoder kegg 2.6 - primed + droput 20%
se.kegg2.6 <- readFeatures("results/TCGA_gexp_norm/kegg_v2.6/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.6) <- as.character(paths[, 1])

brca2.6 <- se.kegg2.6[, se.kegg2.6$project_id == "TCGA-BRCA" & se.kegg2.6$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.6))
paths.lm2.6 <- lmFit(assay(brca2.6), mod)
paths.lm2.6 <- eBayes(paths.lm2.6)
paths.top2.6 <- topTable(paths.lm2.6, n = Inf, coef = 2)
paths.top2.6[cancer.kegg, ]
sum(paths.top2.6[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.6 <- sapply(tumors, computeTPR, se =  se.kegg2.6)

paths2.6 <- lapply(tumors, getTumorPaths, se =  se.kegg2.6)
pathMat2.6 <- sapply(paths2.6, function(x) x$adj.P.Val < 0.05)
colnames(pathMat2.6) <- tumors
rownames(pathMat2.6) <- cancer.keg

## Load autoencoder kegg 2.7 - primed + droput 50%
se.kegg2.7 <- readFeatures("results/TCGA_gexp_norm/kegg_v2.7/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.7) <- as.character(paths[, 1])

brca2.7 <- se.kegg2.7[, se.kegg2.7$project_id == "TCGA-BRCA" & se.kegg2.7$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.7))
paths.lm2.7 <- lmFit(assay(brca2.7), mod)
paths.lm2.7 <- eBayes(paths.lm2.7)
paths.top2.7 <- topTable(paths.lm2.7, n = Inf, coef = 2)
paths.top2.7[cancer.kegg, ]
sum(paths.top2.7[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.7 <- sapply(tumors, computeTPR, se =  se.kegg2.7)
paths2.7 <- lapply(tumors, getTumorPaths, se =  se.kegg2.7)
pathMat2.7 <- sapply(paths2.7, function(x) x$adj.P.Val < 0.05)
colnames(pathMat2.7) <- tumors
rownames(pathMat2.7) <- cancer.keg


## Load autoencoder kegg 2.7b - primed + droput 50%
se.kegg2.7b <- readFeatures("results/TCGA_gexp_norm/kegg_v2.7b/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.7b) <- as.character(paths[, 1])

## Load autoencoder kegg 2.8 - primed + droput 75%
se.kegg2.8 <- readFeatures("results/TCGA_gexp_norm/kegg_v2.8/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.8) <- as.character(paths[, 1])

brca2.8 <- se.kegg2.8[, se.kegg2.8$project_id == "TCGA-BRCA" & se.kegg2.8$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.8))
paths.lm2.8 <- lmFit(assay(brca2.8), mod)
paths.lm2.8 <- eBayes(paths.lm2.8)
paths.top2.8 <- topTable(paths.lm2.8, n = Inf, coef = 2)
paths.top2.8[cancer.kegg, ]
sum(paths.top2.8[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.8 <- sapply(tumors, computeTPR, se =  se.kegg2.8)
paths2.8 <- lapply(tumors, getTumorPaths, se =  se.kegg2.8)
pathMat2.8 <- sapply(paths2.8, function(x) x$adj.P.Val < 0.05)
colnames(pathMat2.8) <- tumors
rownames(pathMat2.8) <- cancer.keg

## Load autoencoder kegg 2.8 - primed + droput 90%
se.kegg2.9 <- readFeatures("results/TCGA_gexp_norm/kegg_v2.9/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.9) <- as.character(paths[, 1])

brca2.9 <- se.kegg2.9[, se.kegg2.9$project_id == "TCGA-BRCA" & se.kegg2.9$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.9))
paths.lm2.9 <- lmFit(assay(brca2.9), mod)
paths.lm2.9 <- eBayes(paths.lm2.9)
paths.top2.9 <- topTable(paths.lm2.9, n = Inf, coef = 2)
paths.top2.9[cancer.kegg, ]
sum(paths.top2.9[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.9 <- sapply(tumors, computeTPR, se =  se.kegg2.9)



## Load autoencoder kegg 3
se.kegg3 <- readFeatures("results/TCGA_gexp_norm/kegg_v3.1/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg3) <- as.character(paths[, 1])

brca3 <- se.kegg3[, se.kegg3$project_id == "TCGA-BRCA" & se.kegg3$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca3))
paths.lm3 <- lmFit(assay(brca3), mod)
paths.lm3 <- eBayes(paths.lm3)
paths.top3 <- topTable(paths.lm3, n = Inf, coef = 2)
paths.top3[cancer.kegg, ]

sum(paths.top3[cancer.kegg, ]$adj.P.Val < 0.05)
tpr3 <- sapply(tumors, computeTPR, se =  se.kegg3)
paths3 <- lapply(tumors, getTumorPaths, se =  se.kegg3)
pathMat3 <- sapply(paths3, function(x) x$adj.P.Val < 0.05)
colnames(pathMat3) <- tumors
rownames(pathMat3) <- cancer.keg

brca3.cont <- brca3[, brca3$sample_type == "Solid Tissue Normal"]
brca3.cont$cot <- sample(c("A", "B"), ncol(brca3.cont), replace = TRUE)

# pc.brca3.cont <- makePCdf(brca3.cont, vars = c("project_id", "sample_type", "cot"))
# png("figures/TCGA_BRCA.cont.kegg.pca.png")
# makePCplot(pc.brca3.cont, "cot")
# dev.off()

mod.fpr <- model.matrix(~cot, colData(brca3.cont))
paths.lm3fpr <- lmFit(assay(brca3.cont ), mod.fpr) %>% eBayes()
paths.top3fpr <- topTable(paths.lm3fpr, n = Inf, coef = 2)
min(paths.top3fpr$P.Value)


vsd.brca3.cont <- vsd.ori[,  vsd.ori$project_id == "TCGA-BRCA" & vsd.ori$sample_type == "Solid Tissue Normal"]
pc.vsd.brca3.cont <- makePCdf(vsd.brca3.cont, vars = c("project_id", "sample_type"))
png("figures/TCGA_BRCA.cont.vsd.pca.png")
makePCplot(pc.vsd.brca3.cont, "sample_type")
dev.off()


png("figures/TCGA.kegg_comb.tpr.png")
data.frame(TPR = c(tpr3, tpr2.6, tpr2.7)/length(cancer.kegg), mode = rep(c("primed", "dropout 20%", "dropout 50%"), each = length(tpr3))) %>%
  mutate(mode = factor(mode, levels = c("primed", "dropout 20%", "dropout 50%") )) %>%
  ggplot(aes(x = mode, y = TPR)) + geom_boxplot() +
  theme_bw()
dev.off()

getNfeatures <- function(se){
  sum(rowSds(assay(se)) > 1e-5)
}

png("figures/TCGA.kegg_comb.Nfeatures.png")
data.frame(N = sapply(list(se.kegg2.6, se.kegg2.7, se.kegg3), getNfeatures), mode = c( "dropout 20%", "dropout 50%", "primed")) %>%

  mutate(mode = factor(mode, levels = c("primed", "dropout 20%", "dropout 50%") )) %>%
  ggplot(aes(x = mode, y = N)) + geom_bar(stat = "identity") +
  theme_bw()
dev.off()


## Correlation with N features
kegg.map <- read.table("results/preprocess/kegg_gene_map.tsv", header = TRUE)
kegg.N <- table(kegg.map$PathwayID)
png("figures/TCGA.kegg_comb.TPRvsNgenes.png", height = 300)

Reduce(rbind, lapply(list(pathMat3, pathMat2.6, pathMat2.7), function(x) {
  data.frame(N_genes = as.numeric(table(kegg.map$PathwayID)[cancer.kegg]), TP = rowSums(x) )
}))  %>%
mutate(model = factor(rep(c("primed", "dropout 20%", "dropout 50%"), each = length(cancer.kegg)),
 levels = c("primed", "dropout 20%", "dropout 50%") ))   %>%
 ggplot(aes(x = N_genes, y = TP)) +
 geom_point() +
 theme_bw() +
 geom_smooth(method = "lm") +
 facet_grid(~model)
dev.off()

cor.drop <- cor(t(assay(se.kegg2.6)), t(assay(se.kegg2.7)), method  = "spearman")
cor.drop2 <- cor(t(assay(se.kegg2.8)), t(assay(se.kegg2.7)), method  = "spearman")

cor.rep <- cor(t(assay(se.kegg2.7)), t(assay(se.kegg2.7b)), method  = "spearman")


png("figures/TCGA.kegg_comb.cor_dropout.png")

data.frame(R = diag(cor.drop), N_genes = as.numeric(kegg.N[rownames(se.kegg2.6)])) %>%
filter(N_genes < 600) %>%
  ggplot(aes(x = N_genes, y = abs(R))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
dev.off()

summary(abs(diag(cor.drop)))

cor.primed <- cor(t(assay(se.kegg3)), t(assay(se.kegg2.7)), method  = "spearman")
cor.primed2 <- cor(t(assay(se.kegg3)), t(assay(se.kegg2.8)), method  = "spearman")


png("figures/TCGA.kegg_comb.cor_primed.png")

data.frame(R = diag(cor.primed), N_genes = as.numeric(kegg.N[rownames(se.kegg2.7)])) %>%
filter(N_genes < 600) %>%
  ggplot(aes(x = N_genes, y = abs(R))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
dev.off()

summary(abs(diag(cor.primed)))

## Load autoencoder hipathia 3
se.hip3 <- readFeatures("results/TCGA_gexp_norm/hipathia_v3.1/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
paths2 <- read.table("results/TCGA_gexp_norm/hipathia_v3.1/model_trained/pathways_names.txt", header = TRUE, sep = "\t")

rownames(se.hip3) <- c(as.character(paths2[, 1]), "NA")


brca3.hip <- se.hip3[, se.hip3$project_id == "TCGA-BRCA" & se.hip3$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca3.hip))
paths.lm3h <- lmFit(assay(brca3.hip ), mod)
paths.lm3h <- eBayes(paths.lm3h)
paths.top3h <- topTable(paths.lm3h, n = Inf, coef = 2)

tpr.hip.ml <- sapply(tumors, computeTPRhipathia, se =  se.hip3[rownames(se.hip3) != "NA",])

## Hipathia
trans_data <- translate_data(vsd.ori, "hsa")
exp_data <- normalize_data(trans_data)
pathways <- load_pathways(species = "hsa")

hip.res <- hipathia(exp_data, pathways, decompose = FALSE, verbose = TRUE)
save(hip.res, file = "results/TCGA_gexp/hipathia.res.Rdata")

## Differential expression
path_vals <- get_paths_data(hip.res)

brca.hip <- path_vals[, path_vals$project_id == "TCGA-BRCA" & path_vals$sample_type != "Metastatic"]

sample_group <- brca3$sample_type
comp <- do_wilcoxon(brca.hip, sample_group, g1 = "Primary Tumor", g2 = "Solid Tissue Normal")
pathways_summary <- get_pathways_summary(comp, pathways)

tpr.hip <- sapply(tumors, computeTPRhipathia, se =  path_vals)

save(tpr.hip, tpr.hip.ml, file =  "results/TCGA_gexp/hipathia.TPR.Rdata")

png("figures/TCGA.kegg_hiapthia.TPR.png")
data.frame(TPR = c(colSums(tpr.hip.ml > 0), colSums(tpr.hip > 0))/14, method = rep(c("DNN kegg", "Hipathia"), each = ncol(tpr.hip.ml ))) %>%
  ggplot(aes(x = method, y = TPR)) + geom_boxplot() +
  theme_bw()
dev.off()

png("figures/TCGA.kegg_hipathia.TPR50.png")
data.frame(TPR = c(colSums(tpr.hip.ml > 50), colSums(tpr.hip > 50))/14, method = rep(c("DNN kegg", "Hipathia"), each = ncol(tpr.hip.ml ))) %>%
  ggplot(aes(x = method, y = TPR)) + geom_boxplot() +
  theme_bw()
dev.off()

png("figures/TCGA.kegg_hipathia.median.png")
data.frame(median = c(colMedians(tpr.hip.ml), colMedians(tpr.hip )), method = rep(c("DNN kegg", "Hipathia"), each = ncol(tpr.hip.ml ))) %>%
  ggplot(aes(x = method, y = median)) + geom_boxplot() +
  theme_bw()
dev.off()




hip.paths <- hip.res[["paths"]][rowSds(assay(hip.res[["paths"]])) > 1e-5,]
hiprows <- unique(gsub("P-", "", rownames(hip.paths)))
hiprows <- paste0("path:", unique(gsub("-.*$", "", hiprows)))

## Correlation with combined pathways
cors <- cor(t(assay(se.kegg3[hiprows, ])), t(assay(hip.paths  )))
heatmap(cors,  col = cm.colors(256))

max.cors <- lapply(rownames(cors), function(x) {
  vec <- cors[x, ]
  p <- gsub("path:", "", x)

  vec.mini <- vec[grep(p, colnames(cors))]
  data.frame(path.cor = vec.mini[which.max(abs(vec.mini))], top.path = names(vec.mini)[which.max(abs(vec.mini))],
    max.cor = vec[which.max(abs(vec))], max.path = colnames(cors)[which.max(abs(vec))])
})
max.cors <- Reduce(rbind, max.cors)
max.cors$ml_path <- rownames(se.kegg3[hiprows, ])
max.cors$n_paths <- sapply(max.cors$ml_path, function(x)  length(grep(gsub("path:", "", x), colnames(cors))))


ggplot(max.cors, aes(x = path.cor, y = max.cor)) + geom_point()
ggplot(max.cors, aes(x = path.cor)) + geom_histogram()



## Correlation with effective pathways
cors.eff <- cor(t(assay(se.hip3)), t(assay(hip.res[["paths"]])))

## Compute correlations
png("figures/TCGA.hip_dnn.corr.png")
hist(diag(cors.eff))
dev.off()


## Create random variables
png("figures/TCGA.hip_dnn.random.corr.png")
hist(c(cors.eff[upper.tri(cors.eff)], cors.eff[lower.tri(cors.eff)]))
dev.off()

cancer <- unique(se.hip3$project_id)
cors.eff.cancer <- lapply(cancer, function(x){
  cor(t(assay(se.hip3[,se.hip3$project_id == x])), t(assay(path_vals[, path_vals$project_id == x])))
})

png("figures/TCGA.hip_dnn.project.corr.png")

data.frame(cors = unlist(lapply(cors.eff.cancer, diag)), project = rep(cancer, each = ncol(cors.eff ))) %>%
  ggplot(aes(x = cors)) +
  facet_wrap(~project) +
  geom_histogram() +
  theme_bw()

dev.off()

# png("pac.png")
# data.frame(ml = as.numeric(assay(se.hip3[27,])),
#             hip = as.numeric(assay(hip.res[["paths"]][27, ])),
#             cancer = se.hip3$project_id) %>%
#             ggplot(aes( x = ml, y = hip, color = cancer) ) +
#             geom_point()
# dev.off()

cancer <- unique(se.hip3$project_id)
sapply(cancer, function(x){
  cor(as.numeric(assay(se.hip3[27,se.hip3$project_id == x])), as.numeric(assay(hip.res[["paths"]][27, se.hip3$project_id == x])))

})

sapply(cancer, function(x){
  cor(as.numeric(assay(se.hip3[27,se.hip3$project_id == x])), as.numeric(assay(hip.res[["paths"]][2, se.hip3$project_id == x])))

})

cor.df <- data.frame(cor = diag(cors.eff),
          path = gsub("-.*$", "", gsub("P-", "", rownames(cors.eff))))

png("pac.png", width = 2000, height = 2000)
 cor.df %>%
            ggplot(aes( x = cor) ) +
            geom_histogram() +
            facet_wrap( ~ path, scales = "free_y") +
            theme_bw()
dev.off()


## Correlation with autoencoder
se.auto <- readFeatures("results/TCGA_gexp_norm/autoencod_v2.1/model_features/dense.tsv", vsd.ori)

cor.auto <- cor(t(assay(se.kegg3)), t(assay(se.auto  )))



## Initial version: models without primed initialization
## Load autoencoder kegg 1
se.kegg1 <- readFeatures("results/TCGA_gexp_norm/kegg_v1.3/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg1) <- as.character(paths[, 1])

brca1 <- se.kegg1[, se.kegg1$project_id == "TCGA-BRCA" & se.kegg1$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca1))
paths.lm1 <- lmFit(assay(brca1), mod)
paths.lm1 <- eBayes(paths.lm1)
paths.top1 <- topTable(paths.lm1, n = Inf, coef = 2)
paths.top1[cancer.kegg, ]
sum(paths.top1[cancer.kegg, ]$adj.P.Val < 0.05)

makeHeatmap(brca1[cancer.kegg, ])

tpr1 <- sapply(tumors, computeTPR, se =  se.kegg1)



## Load autoencoder kegg 2.1
se.kegg2.1 <- readFeatures("results/TCGA_gexp_kegg/v2.1/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.1) <- as.character(paths[, 1])

brca2.1 <- se.kegg2.1[, se.kegg2.1$project_id == "TCGA-BRCA" & se.kegg2.1$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.1))
paths.lm2.1 <- lmFit(assay(brca2.1), mod)
paths.lm2.1 <- eBayes(paths.lm2.1)
paths.top2.1 <- topTable(paths.lm2.1, n = Inf, coef = 2)
paths.top2.1[cancer.kegg, ]
sum(paths.top2.1[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.1 <- sapply(tumors, computeTPR, se =  se.kegg2.1)

## Load autoencoder kegg 2.2
se.kegg2.2 <- readFeatures("results/TCGA_gexp_kegg/v2.2/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.2) <- as.character(paths[, 1])

brca2.2 <- se.kegg2.2[, se.kegg2.2$project_id == "TCGA-BRCA" & se.kegg2.2$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.2))
paths.lm2.2 <- lmFit(assay(brca2.2), mod)
paths.lm2.2 <- eBayes(paths.lm2.2)
paths.top2.2 <- topTable(paths.lm2.2, n = Inf, coef = 2)
paths.top2.2[cancer.kegg, ]
sum(paths.top2.2[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.2 <- sapply(tumors, computeTPR, se =  se.kegg2.2)



## Load autoencoder kegg 2.4
se.kegg2.4 <- readFeatures("results/TCGA_gexp_kegg/v2.4/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
rownames(se.kegg2.4) <- as.character(paths[, 1])

brca2.4 <- se.kegg2.4[, se.kegg2.4$project_id == "TCGA-BRCA" & se.kegg2.4$sample_type != "Metastatic"]

mod <- model.matrix(~sample_type, colData(brca2.4))
paths.lm2.4 <- lmFit(assay(brca2.4), mod)
paths.lm2.4 <- eBayes(paths.lm2.4)
paths.top2.4 <- topTable(paths.lm2.4, n = Inf, coef = 2)
paths.top2.4[cancer.kegg, ]
sum(paths.top2.4[cancer.kegg, ]$adj.P.Val < 0.05)
tpr2.4 <- sapply(tumors, computeTPR, se =  se.kegg2.4)
