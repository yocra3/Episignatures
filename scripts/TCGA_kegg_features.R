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

  sub.se <- se[, se$project_id == tumor & se.kegg1$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  mod <- model.matrix(~sample_type, colData(sub.se))
  paths.lm <- lmFit(assay(sub.se), mod)
  paths.lm <- eBayes(paths.lm)
  paths.top <- topTable(paths.lm, n = Inf, coef = 2)
  sum(paths.top[cancer.kegg, ]$adj.P.Val < 0.05)


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
auto.vals <- read_training("results/TCGA_gexp_norm/autoencod_v2.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "autoencoder")
raw.vals <- read_training("results/TCGA_gexp_norm/kegg_v1.3/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - base")
drop20.vals <- read_training("results/TCGA_gexp_norm/kegg_v2.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - dropout 20%")
drop50.vals <- read_training("results/TCGA_gexp_norm/kegg_v2.2/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - dropout 50%")
primed.vals <- read_training("results/TCGA_gexp_norm/kegg_v3.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "kegg - primed")
hipathia.vals <- read_training("results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm_training_evaluation.txt", "hipathia - primed")

all.vals <- Reduce(rbind, list(auto.vals, raw.vals, drop20.vals, drop50.vals, primed.vals,hipathia.vals))



png("figures/TCGA.kegg.training.png", width = 1200)
ggplot(all.vals, aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ model) +
  theme_bw()
dev.off()

png("figures/TCGA.kegg.training2.png", width = 1200)
filter(all.vals, mse < 0.005) %>%
ggplot(aes(x = epoch, y = mse, color = measure, group = measure)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ model) +
  theme_bw()
dev.off()


## Load vsd data ####
vsd.ori <- loadHDF5SummarizedExperiment("results/TCGA_gexp/", prefix = "vsd_norm")

cancer.keg <- c("hsa04010", "hsa04310", "hsa04350", "hsa04370", "hsa04630", "hsa04024", "hsa04151", "hsa04150", "hsa04110", "hsa04210", "hsa04115", "hsa04510", "hsa04520", "hsa03320")
cancer.kegg <- paste0("path:", cancer.keg)
paths <- read.table("results/TCGA_gexp_kegg/v1.1/model_trained/pathways_names.txt", header = TRUE)

tumors <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-THCA", "TCGA-UCEC")


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
data.frame(TPR = c(tpr1, tpr2.1, tpr2.2, tpr3)/length(cancer.kegg), mode = rep(c("base", "dropout 20%", "dropout 50%", "primed"), each = length(tpr1))) %>%
  mutate(mode = factor(mode, levels = c("base", "dropout 20%", "dropout 50%", "primed") )) %>%
  ggplot(aes(x = mode, y = TPR)) + geom_boxplot() +
  theme_bw()
dev.off()

getNfeatures <- function(se){
  sum(rowSds(assay(se)) > 1e-5)
}

png("figures/TCGA.kegg_comb.Nfeatures.png")
data.frame(N = sapply(list(se.kegg1, se.kegg2.1, se.kegg2.2, se.kegg3), getNfeatures), mode = c("base", "dropout 20%", "dropout 50%", "primed")) %>%
  mutate(mode = factor(mode, levels = c("base", "dropout 20%", "dropout 50%", "primed") )) %>%
  ggplot(aes(x = mode, y = N)) + geom_bar(stat = "identity") +
  theme_bw()
dev.off()

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
