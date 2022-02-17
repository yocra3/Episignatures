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


makePCdf <- function(seobj){
  pc <- prcomp(t(assay(seobj)))
  pcdf <- data.frame(pc$x[, 1:10])
  pcdf <- cbind(pcdf, colData(seobj)[, c("diagnosis2", "age", "sex")])
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
    diagnosis2 = c("Control" = "black", "cCD" = "darkgreen", "UC" = "blue",
                "iCD" = "lightgreen"),
    sex = c("Female" = "purple", "Male" = "lightblue")
  )

  pheatmap(assay(seObj), scale = "row",
           annotation_col  = data.frame(colData(seObj)[, c("diagnosis2", "sex"), drop = FALSE]),
           annotation_colors =  col_colors,
           show_rownames = FALSE)

}
trainSVM <-  function(seobj){
  mat <- data.matrix(assay(seobj))
  df <- data.frame(pathClass = factor(seobj$diagnosis2), t(mat))
  model_svm <- svm(pathClass ~ ., df)
}

getPathwayCor <- function(path){
  path2 <- gsub(".", ":", path, fixed = TRUE)
  genes <- subset(path_map, PathwayID == path2)$Symbol
  pc <- prcomp(t(assay(vsd.uc[rownames(vsd.uc) %in% genes,])), rank. = 1)$x[, 1]
  cor(pc, t(assay(se.kegg.uc[path,])))

}





## Load vsd data ####
vsd.ori <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm_TCGAgenes_")
vsd.ori$diagnosis2 <- relevel(factor(vsd.ori$diagnosis2), ref = "Control")
## Create train and test indexes
sample <- sample.split(vsd.ori$diagnosis2, SplitRatio = .80)


## Get features associated with CD
mod <- model.matrix(~diagnosis2 + age + sex, colData(vsd.ori))

lmori <- lmFit(assay(vsd.ori), mod) %>% eBayes()
tab.ori <- topTable(lmori, coef = 2:4, n = Inf)
featsori <- rownames(subset(tab.ori, adj.P.Val  < 0.05))

## See age
tab.ori.age <- topTable(lmori, coef = 5, n = Inf)
featsori.age <- rownames(subset(tab.ori.age, adj.P.Val  < 0.05))


## See sex
tab.ori.sex <- topTable(lmori, coef = 6, n = Inf)
featsori.sex <- rownames(subset(tab.ori.sex, adj.P.Val  < 0.05))

## Make PCs
pc.ori <- makePCdf(vsd.ori)
pc.ori.feats <- makePCdf(vsd.ori[featsori, ])

png("figures/SRP042228.raw.pca.png")
makePCplot(pc.ori, "diagnosis2")
dev.off()

makePCplot(pc.ori, "age")
makePCplot(pc.ori, "sex")

makePCplot(pc.ori.feats, "diagnosis2")
makePCplot(pc.ori.feats, "age")
makePCplot(pc.ori.feats, "sex")

## Heatmaps
#makeHeatmap(vsd.ori)
makeHeatmap(vsd.ori[featsori, ])

## SVM
svm.ori <- trainSVM(vsd.ori[featsori, sample])
pred.ori <- predict(svm.ori, t(assay(vsd.ori[featsori, !sample])))
table(prediction = pred.ori, real = vsd.ori$diagnosis2[!sample] )

## UC
vsd.uc <- vsd.ori[, vsd.ori$diagnosis != "CD"]
vsd.uc$diagnosis2 <- droplevels(vsd.uc$diagnosis2)
mod.uc <- model.matrix(~diagnosis2 + age + sex, colData(vsd.uc ))
lm.uc <- lmFit(assay(vsd.uc), mod.uc) %>% eBayes()
tab.uc <- topTable(lm.uc, coef = 2, n = Inf)
feats.uc  <- rownames(subset(tab.uc  , adj.P.Val  < 0.1))

## DNN kegg
se.kegg <- readFeatures("results/SRP042228/kegg_v3.1/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
paths <- read.table("results/TCGA_gexp_norm/kegg_v3.1/model_trained/pathways_names.txt", header = TRUE, sep = "\t")
rownames(se.kegg)  <- gsub(":", ".",as.character( paths[, 1]))

## SVM
svm.kegg <- trainSVM(se.kegg[, sample])
pred.kegg <- predict(svm.kegg, t(assay(se.kegg[, !sample])))
table(prediction = pred.kegg, real = vsd.ori$diagnosis2[!sample] )

lm.kegg <- lmFit(assay(se.kegg), mod) %>% eBayes()
tab.kegg <- topTable(lm.kegg, coef = 2:4, n = Inf)
feats.kegg  <- rownames(subset(tab.kegg , adj.P.Val  < 0.05))

svm.kegg.f <- trainSVM(se.kegg[feats.kegg, sample])
pred.kegg.f <- predict(svm.kegg.f, t(assay(se.kegg[feats.kegg, !sample])))
table(prediction = pred.kegg.f, real = vsd.ori$diagnosis2[!sample] )


pc.kegg <- makePCdf(se.kegg)

png("figures/SRP042228.kegg_all_dnn.pca.png")
makePCplot(pc.kegg, "diagnosis2")
dev.off()


se.kegg.uc <- se.kegg[, se.kegg$diagnosis != "CD"]
se.kegg.uc$diagnosis2 <- droplevels(se.kegg.uc$diagnosis2)
mod.uc <- model.matrix(~diagnosis2 + age + sex, colData(se.kegg.uc ))
lm.kegg.uc <- lmFit(assay(se.kegg.uc), mod.uc) %>% eBayes()
tab.kegg.uc <- topTable(lm.kegg.uc , coef = 2, n = Inf)
feats.kegg.uc  <- rownames(subset(tab.kegg.uc  , adj.P.Val  < 0.05))


png("figures/SRP042228.kegg_uc.heatmap.png")
makeHeatmap(se.kegg.uc[feats.kegg.uc , ])
dev.off()

png("figures/SRP042228.kegg_uc.heatmap_all.png")
makeHeatmap(se.kegg.uc)
dev.off()

path_map <- read.table("results/preprocess/kegg_gene_map.tsv", header = TRUE )
genes1 <- subset(path_map, PathwayID == "path:hsa01521")$Symbol


png("figures/SRP042228.vsd_uc.kegg_hsa0152.heatmap.png")
makeHeatmap(vsd.uc[genes1 , ])
dev.off()

cot <- prcomp(t(assay(vsd.uc[genes1,])))
png("figures/SRP042228.vsd_uc.kegg_hsa0152.cor.png")
plot(cot$x[, 1], t(assay(se.kegg.uc["path.hsa01521",])), col = se.kegg.uc$diagnosis2)
dev.off()
cor(cot$x[, 1], t(assay(se.kegg.uc["path.hsa01521",])))
summary(lm(cot$x[, 1] ~ vsd.uc$diagnosis2 + vsd.uc$sex + vsd.uc$age))
summary(lm(t(assay(se.kegg.uc["path.hsa01521",])) ~ vsd.uc$diagnosis2 + vsd.uc$sex + vsd.uc$age))


genes2 <- subset(path_map, PathwayID == "path:hsa04657")$Symbol
png("figures/SRP042228.vsd_uc.kegg_hsa04657.heatmap.png")
makeHeatmap(vsd.uc[genes2, ] )
dev.off()

cot2 <- prcomp(t(assay(vsd.uc[genes2,])))
png("figures/SRP042228.vsd_uc.kegg_hsa04657.cor.png")
plot(cot2$x[, 1], t(assay(se.kegg.uc["path.hsa04657",])), col = se.kegg.uc$diagnosis2)
dev.off()

cor(cot2$x[, 1], t(assay(se.kegg.uc["path.hsa04657",])))

summary(lm(cot2$x[, 1] ~ vsd.uc$diagnosis2 + vsd.uc$sex + vsd.uc$age))
summary(lm(t(assay(se.kegg.uc["path.hsa04657",])) ~ vsd.uc$diagnosis2 + vsd.uc$sex + vsd.uc$age))

path.cores <- sapply(rownames(se.kegg.uc), getPathwayCor)
names(path.cores) <- rownames(se.kegg.uc)
path_n <- sapply(rownames(se.kegg.uc), function(x) sum(gsub(":", ".", path_map$PathwayID) == x))
cor(path_n, abs(path.cores ))

png("figures/SRP042228.uc.vsd_genes.kegg.cor.png")
hist(path.cores, breaks = 40)
dev.off()

path_n2 <- path_n
path_n2[path_n2 > 500] <- 500
png("figures/SRP042228.uc.vsd_genes.kegg.cor_vs_n.png")
plot(path.cores, path_n2)
dev.off()


## DNN hipathia
se.kegg.hip <- readFeatures("results/SRP042228/hipathia_v3.1/model_features/prune_low_magnitude_dense.tsv", vsd.ori)
paths2 <- read.table("results/TCGA_gexp_norm/hipathia_v3.1/model_trained/pathways_names.txt", header = TRUE, sep = "\t")

rownames(se.kegg.hip) <- c(as.character(paths2[, 1]), "NA")

## SVM
svm.kegg <- trainSVM(se.kegg.hip[-1876, sample])
se.kegg.hip2 <- se.kegg.hip
rownames(se.kegg.hip2) <- gsub("-", ".", rownames(se.kegg.hip2))
rownames(se.kegg.hip2) <- gsub(" ", ".", rownames(se.kegg.hip2))

lm.kegg <- lmFit(assay(se.kegg.hip), mod) %>% eBayes()
tab.kegg <- topTable(lm.kegg, coef = 2:4, n = Inf)
feats.kegg  <- rownames(subset(tab.kegg , adj.P.Val  < 0.05))

tab.kegg <- topTable(lm.kegg, coef = 2:4, n = Inf)


pred.kegg <- predict(svm.kegg, t(assay(se.kegg.hip2[-1876, !sample])))
table(prediction = pred.kegg, real = vsd.ori$diagnosis2[!sample] )

svm.kegg.f <- trainSVM(se.kegg.hip[feats.kegg, sample])
feats.kegg2 <- gsub(" ", ".", gsub("-", ".", feats.kegg))
pred.kegg.f <- predict(svm.kegg.f, t(assay(se.kegg.hip2[feats.kegg2, !sample])))
table(prediction = pred.kegg.f, real = vsd.ori$diagnosis2[!sample] )


pc.kegg <- makePCdf(se.kegg.hip)

png("figures/SRP042228.kegg_dnn.pca.png")
makePCplot(pc.kegg, "diagnosis2")
dev.off()

## Hipathia
trans_data <- translate_data(vsd.ori, "hsa")
exp_data <- normalize_data(trans_data)
pathways <- load_pathways(species = "hsa")

hip.res <- hipathia(exp_data, pathways, decompose = FALSE, verbose = TRUE)
save(hip.res, file = "results/SRP042228/hipathia.res.Rdata")

hip.paths <- get_paths_data(hip.res)

pc.hip <- makePCdf(hip.paths)

png("figures/SRP042228.hipathia.pca.png")
makePCplot(pc.hip, "diagnosis2")
dev.off()


## SVM
lm.hip <- lmFit(assay(hip.paths ), mod) %>% eBayes()
tab.hip <- topTable(lm.hip, coef = 2:4, n = Inf)
feats.hip  <- rownames(subset(tab.hip , adj.P.Val  < 0.05))


svm.hip <- trainSVM(hip.paths[, sample])
hip.paths2 <- hip.paths
rownames(hip.paths2) <- gsub("-", ".", rownames(hip.paths2))
rownames(hip.paths2) <- gsub(" ", ".", rownames(hip.paths2))

pred.hip <- predict(svm.hip, t(assay(hip.paths2[, !sample])))
table(prediction = pred.hip, real = vsd.ori$diagnosis2[!sample] )

svm.hip.f <- trainSVM(hip.paths[feats.hip, sample])
feats.hip2 <- gsub(" ", ".", gsub("-", ".", feats.hip))
pred.hip.f <- predict(svm.hip.f, t(assay(hip.paths2[feats.hip2, !sample])))
table(prediction = pred.hip.f, real = vsd.ori$diagnosis2[!sample] )


## Compute correlations
cors <- cor(t(assay(se.kegg.hip)), t(assay(hip.paths)))
png("figures/SRP042228.hip_dnn.corr.png")
hist(diag(cors))
dev.off()


## Create random variables
png("figures/SRP042228.hip_dnn.random.corr.png")
hist(c(cors[upper.tri(cors)], cors[lower.tri(cors)]))
dev.off()

mean(abs(diag(cors)) > 0.5, na.rm = TRUE)
# [1] 0.3008
mean(abs(c(cors[upper.tri(cors)], cors[lower.tri(cors)])) > 0.5, na.rm = TRUE)
# [1] 0.1196877

mean(sapply(seq_len(ncol(cors)), function(i) which.max(abs(cors[i, ])) == i), na.rm = TRUE)
# [1] 0.06982942

cors.kegg <- cor(t(assay(se.kegg.hip)))
png("figures/SRP042228.dnn.random.corr.png")
hist(cors.kegg[upper.tri(cors.kegg)])
dev.off()

png("figures/SRP042228.dnn.random.corr.heat.png", width = 2000, height = 2000)
pheatmap(cors.kegg, cluster_rows = FALSE,  cluster_cols = FALSE)
dev.off()


### Viejo a partir de aquí
## Load model 1.4 layer1
se.mod1.4_l1 <- readFeatures("results/SRP042228/v1.4.prun.2/model_features/prune_low_magnitude_dense.tsv", vsd.ori)

## Remove features with 0 in all samples
se.mod1.4_l1f <- se.mod1.4_l1[rowSds(assay(se.mod1.4_l1)) >  1e-5, ]
lm.mod1.4_l1 <- lmFit(assay(se.mod1.4_l1f), mod) %>% eBayes()
tab.mod1.4_l1 <- topTable(lm.mod1.4_l1, coef = 2:4, n = Inf)
feats.mod1.4_l1 <- rownames(subset(tab.mod1.4_l1, adj.P.Val  < 0.05))


## Make PCs
pc.mod1.4_l1 <- makePCdf(se.mod1.4_l1f)
pc.mod1.4_l1.feats <- makePCdf(se.mod1.4_l1f[feats.mod1.4_l1, ])

makePCplot(pc.mod1.4_l1, "diagnosis2")
makePCplot(pc.mod1.4_l1, "age")
makePCplot(pc.mod1.4_l1, "sex")

makePCplot(pc.mod1.4_l1.feats, "diagnosis2")
makePCplot(pc.mod1.4_l1.feats, "age")
makePCplot(pc.mod1.4_l1.feats, "sex")

makeHeatmap(se.mod1.4_l1f)
makeHeatmap(se.mod1.4_l1f[feats.mod1.4_l1, ])

## SVM
svm.mod1.4_l1 <- trainSVM(se.mod1.4_l1f[, sample])
pred.mod1.4_l1 <- predict(svm.mod1.4_l1, t(assay(se.mod1.4_l1f[, !sample])))
table(prediction = pred.mod1.4_l1, real = vsd.ori$diagnosis2[!sample] )

svm.mod1.4_l1f <- trainSVM(se.mod1.4_l1f[feats.mod1.4_l1, sample])
pred.mod1.4_l1f <- predict(svm.mod1.4_l1f, t(assay(se.mod1.4_l1f[feats.mod1.4_l1, !sample])))
table(prediction = pred.mod1.4_l1f, real = vsd.ori$diagnosis2[!sample] )


## Load model sex 1.4 layer1
se.mod1.4sex_l1 <- readFeatures("results/SRP042228/v1.6.prun.2.sex/model_features/prune_low_magnitude_dense.tsv", vsd.ori)

## Remove features with 0 in all samples
se.mod1.4sex_l1f <- se.mod1.4sex_l1[rowSds(assay(se.mod1.4sex_l1)) > 1e-5, ]
lm.mod1.4sex_l1 <- lmFit(assay(se.mod1.4sex_l1f), mod) %>% eBayes()
tab.mod1.4sex_l1 <- topTable(lm.mod1.4sex_l1, coef = 2:4, n = Inf)
feats.mod1.4sex_l1 <- rownames(subset(tab.mod1.4sex_l1, adj.P.Val  < 0.05))

## See age
tab.mod1.4sex_l1.age <- topTable(lm.mod1.4sex_l1, coef = 5, n = Inf)
feats.mod1.4sex_l1.age <- rownames(subset(tab.mod1.4sex_l1.age, adj.P.Val  < 0.05))

## See sex
tab.mod1.4sex_l1.sex <- topTable(lm.mod1.4sex_l1, coef = 6, n = Inf)
feats.mod1.4sex_l1.sex <- rownames(subset(tab.mod1.4sex_l1.sex, adj.P.Val  < 0.05))

makeHeatmap(se.mod1.4sex_l1f)
makeHeatmap(se.mod1.4sex_l1f[feats.mod1.4sex_l1, ])
makeHeatmap(se.mod1.4sex_l1f[feats.mod1.4sex_l1.sex, ])




## control vs UC
tab.mod1.4sex_l1.uc <- topTable(lm.mod1.4sex_l1, coef = 4, n = Inf)


## SVM
svm.mod1.4sex_l1 <- trainSVM(se.mod1.4sex_l1f[, sample])
pred.mod1.4sex_l1 <- predict(svm.mod1.4sex_l1, t(assay(se.mod1.4sex_l1f[, !sample])))
table(prediction = pred.mod1.4sex_l1, real = vsd.ori$diagnosis2[!sample] )

svm.mod1.4sex_l1f <- trainSVM(se.mod1.4sex_l1f[feats.mod1.4sex_l1, sample])
pred.mod1.4sex_l1f <- predict(svm.mod1.4sex_l1f, t(assay(se.mod1.4sex_l1f[feats.mod1.4sex_l1, !sample])))
table(prediction = pred.mod1.4sex_l1f, real = vsd.ori$diagnosis2[!sample] )

## SVM for CD
svm.mod1.4sex_l1cd <- trainSVM(se.mod1.4sex_l1f[, sample & se.mod1.4sex_l1f$diagnosis == "CD"])
pred.mod1.4sex_l1cd <- predict(svm.mod1.4sex_l1cd, t(assay(se.mod1.4sex_l1f[, !sample  & se.mod1.4sex_l1f$diagnosis == "CD"])))
table(prediction = pred.mod1.4sex_l1cd, real = vsd.ori$diagnosis2[!sample  & se.mod1.4sex_l1f$diagnosis == "CD"] )

tab.mod1.4sex_cd <- topTable(lm.mod1.4sex_l1, coef = 3, n = Inf)
feats.mod1.4sex_cd <- rownames(subset(tab.mod1.4sex_cd, adj.P.Val  < 0.05))

svm.mod1.4sex_l1fcd <- trainSVM(se.mod1.4sex_l1f[feats.mod1.4sex_cd, sample & se.mod1.4sex_l1f$diagnosis == "CD"])
pred.mod1.4sex_l1fcd <- predict(svm.mod1.4sex_l1fcd, t(assay(se.mod1.4sex_l1f[feats.mod1.4sex_cd, !sample & se.mod1.4sex_l1f$diagnosis == "CD"])))
table(prediction = pred.mod1.4sex_l1fcd, real = vsd.ori$diagnosis2[!sample & se.mod1.4sex_l1f$diagnosis == "CD"] )

## SVM for control
svm.mod1.4sex_l1cn <- trainSVM(se.mod1.4sex_l1f[, sample & se.mod1.4sex_l1f$diagnosis != "CD"])
pred.mod1.4sex_l1cn <- predict(svm.mod1.4sex_l1cn, t(assay(se.mod1.4sex_l1f[, !sample  & se.mod1.4sex_l1f$diagnosis != "CD"])))
table(prediction = pred.mod1.4sex_l1cn, real = vsd.ori$diagnosis2[!sample  & se.mod1.4sex_l1f$diagnosis != "CD"] )

svm.mod1.4sex_l1fcd <- trainSVM(se.mod1.4sex_l1f[feats.mod1.4sex_l1, sample & se.mod1.4sex_l1f$diagnosis == "CD"])
pred.mod1.4sex_l1fcd <- predict(svm.mod1.4sex_l1fcd, t(assay(se.mod1.4sex_l1f[feats.mod1.4sex_l1, !sample & se.mod1.4sex_l1f$diagnosis == "CD"])))
table(prediction = pred.mod1.4sex_l1fcd, real = vsd.ori$diagnosis2[!sample & se.mod1.4sex_l1f$diagnosis == "CD"] )

## Load model 1.4 layer2
se.mod1.4_l2 <- readFeatures("results/SRP042228/v1.4/model_features/dense_1.tsv", vsd.ori)

## Remove features with 0 in all samples
se.mod1.4_l2f <- se.mod1.4_l2[rowSds(assay(se.mod1.4_l2)) > 0, ]
lm.mod1.4_l2 <- lmFit(assay(se.mod1.4_l2f), mod) %>% eBayes()
tab.mod1.4_l2 <- topTable(lm.mod1.4_l2, coef = 2:4, n = Inf)
feats.mod1.4_l2 <- rownames(subset(tab.mod1.4_l2, adj.P.Val  < 0.05))

## Make PCs
pc.mod1.4_l2 <- makePCdf(se.mod1.4_l2f)
pc.mod1.4_l2.feats <- makePCdf(se.mod1.4_l2f[feats.mod1.4_l2, ])

makePCplot(pc.mod1.4_l2, "diagnosis2")
makePCplot(pc.mod1.4_l2, "age")
makePCplot(pc.mod1.4_l2, "sex")

makePCplot(pc.mod1.4_l2.feats, "diagnosis2")
makePCplot(pc.mod1.4_l2.feats, "age")
makePCplot(pc.mod1.4_l2.feats, "sex")

makeHeatmap(se.mod1.4_l2f)
makeHeatmap(se.mod1.4_l2f[feats.mod1.4_l2, ])

## SVM
svm.mod1.4_l2 <- trainSVM(se.mod1.4_l2f[, sample])
pred.mod1.4_l2 <- predict(svm.mod1.4_l2, t(assay(se.mod1.4_l2f[, !sample])))
table(prediction = pred.mod1.4_l2, real = vsd.ori$diagnosis2[!sample] )

svm.mod1.4_l2f <- trainSVM(se.mod1.4_l2f[feats.mod1.4_l2, sample])
pred.mod1.4_l2f <- predict(svm.mod1.4_l2f, t(assay(se.mod1.4_l2f[feats.mod1.4_l2, !sample])))
table(prediction = pred.mod1.4_l2f, real = vsd.ori$diagnosis2[!sample] )


## Load autoencoder 1.4 layer1
se.auto1.4_l1 <- readFeatures("results/SRP042228/autoencod_autosom_v1.4.prun.2/model_features/prune_low_magnitude_dense.tsv", vsd.ori)

## Remove features with 0 in all samples
se.auto1.4_l1f <- se.auto1.4_l1[rowSds(assay(se.auto1.4_l1)) > 1e-5, ]
lm.auto1.4_l1 <- lmFit(assay(se.auto1.4_l1f), mod) %>% eBayes()
tab.auto1.4_l1 <- topTable(lm.auto1.4_l1, coef = 2:4, n = Inf)
feats.auto1.4_l1 <- rownames(subset(tab.auto1.4_l1, adj.P.Val  < 0.05))

## See age
tab.auto1.4_l1.age <- topTable(lm.auto1.4_l1, coef = 5, n = Inf)
feats.auto1.4_l1.age <- rownames(subset(tab.auto1.4_l1.age, adj.P.Val  < 0.05))


## See sex
tab.auto1.4_l1.sex <- topTable(lm.auto1.4_l1, coef = 6, n = Inf)
feats.auto1.4_l1.sex <- rownames(subset(tab.auto1.4_l1.sex, adj.P.Val  < 0.05))


## Make PCs
pc.auto1.4_l1 <- makePCdf(se.auto1.4_l1f)
pc.auto1.4_l1.feats <- makePCdf(se.auto1.4_l1f[feats.auto1.4_l1, ])

makePCplot(pc.auto1.4_l1, "diagnosis2")
makePCplot(pc.auto1.4_l1, "age")
makePCplot(pc.auto1.4_l1, "sex")

makePCplot(pc.auto1.4_l1.feats, "diagnosis2")
makePCplot(pc.auto1.4_l1.feats, "age")
makePCplot(pc.auto1.4_l1.feats, "sex")

makeHeatmap(se.auto1.4_l1f)
makeHeatmap(se.auto1.4_l1f[feats.auto1.4_l1, ])

## SVM
svm.auto1.4_l1 <- trainSVM(se.auto1.4_l1f[, sample])
pred.auto1.4_l1 <- predict(svm.auto1.4_l1, t(assay(se.auto1.4_l1f[, !sample])))
table(prediction = pred.auto1.4_l1, real = vsd.ori$diagnosis2[!sample] )

svm.auto1.4_l1f <- trainSVM(se.auto1.4_l1f[feats.auto1.4_l1, sample])
pred.auto1.4_l1f <- predict(svm.auto1.4_l1f, t(assay(se.auto1.4_l1f[feats.auto1.4_l1, !sample])))
table(prediction = pred.auto1.4_l1f, real = vsd.ori$diagnosis2[!sample] )


## Load autoencoder 1.4 layer2
se.auto1.4_l2 <- readFeatures("results/SRP042228/autoencod_v1.4/model_features/dense_1.tsv", vsd.ori)

## Remove features with 0 in all samples
se.auto1.4_l2f <- se.auto1.4_l2[rowMeans(assay(se.auto1.4_l2) == 0) < 1, ]
lm.auto1.4_l2 <- lmFit(assay(se.auto1.4_l2f), mod) %>% eBayes()
tab.auto1.4_l2 <- topTable(lm.auto1.4_l2, coef = 2:4, n = Inf)
feats.auto1.4_l2 <- rownames(subset(tab.auto1.4_l2, adj.P.Val  < 0.05))

## Make PCs
pc.auto1.4_l2 <- makePCdf(se.auto1.4_l2f)
pc.auto1.4_l2.feats <- makePCdf(se.auto1.4_l2f[feats.auto1.4_l2, ])

makePCplot(pc.auto1.4_l2, "diagnosis2")
makePCplot(pc.auto1.4_l2, "age")
makePCplot(pc.auto1.4_l2, "sex")

makePCplot(pc.auto1.4_l2.feats, "diagnosis2")
makePCplot(pc.auto1.4_l2.feats, "age")
makePCplot(pc.auto1.4_l2.feats, "sex")

makeHeatmap(se.auto1.4_l2f)
makeHeatmap(se.auto1.4_l2f[feats.auto1.4_l2, ])

## SVM
svm.auto1.4_l2 <- trainSVM(se.auto1.4_l2f[, sample])
pred.auto1.4_l2 <- predict(svm.auto1.4_l2, t(assay(se.auto1.4_l2f[, !sample])))
table(prediction = pred.auto1.4_l2, real = vsd.ori$diagnosis2[!sample] )

svm.auto1.4_l2f <- trainSVM(se.auto1.4_l2f[feats.auto1.4_l2, sample])
pred.auto1.4_l2f <- predict(svm.auto1.4_l2f, t(assay(se.auto1.4_l2f[feats.auto1.4_l2, !sample])))
table(prediction = pred.auto1.4_l2f, real = vsd.ori$diagnosis2[!sample] )


## Load model 1.5 layer1
se.mod1.5_l1 <- readFeatures("results/SRP042228/v1.5/model_features/dense.tsv", vsd.ori)

## Remove features with 0 in all samples
se.mod1.5_l1f <- se.mod1.5_l1[rowMeans(assay(se.mod1.5_l1) == 0) < 1, ]
lm.mod1.5_l1 <- lmFit(assay(se.mod1.5_l1f), mod) %>% eBayes()
tab.mod1.5_l1 <- topTable(lm.mod1.5_l1, coef = 2:4, n = Inf)
feats.mod1.5_l1 <- rownames(subset(tab.mod1.5_l1, adj.P.Val  < 0.05))

## Make PCs
pc.mod1.5_l1 <- makePCdf(se.mod1.5_l1f)
pc.mod1.5_l1.feats <- makePCdf(se.mod1.5_l1f[feats.mod1.5_l1, ])

makePCplot(pc.mod1.5_l1, "diagnosis2")
makePCplot(pc.mod1.5_l1, "age")
makePCplot(pc.mod1.5_l1, "sex")

makePCplot(pc.mod1.5_l1.feats, "diagnosis2")
makePCplot(pc.mod1.5_l1.feats, "age")
makePCplot(pc.mod1.5_l1.feats, "sex")

makeHeatmap(se.mod1.5_l1f)
makeHeatmap(se.mod1.5_l1f[feats.mod1.5_l1, ])

## SVM
svm.mod1.5_l1 <- trainSVM(se.mod1.5_l1f[, sample])
pred.mod1.5_l1 <- predict(svm.mod1.5_l1, t(assay(se.mod1.5_l1f[, !sample])))
table(prediction = pred.mod1.5_l1, real = vsd.ori$diagnosis2[!sample] )

svm.mod1.5_l1f <- trainSVM(se.mod1.5_l1f[feats.mod1.5_l1, sample])
pred.mod1.5_l1f <- predict(svm.mod1.5_l1f, t(assay(se.mod1.5_l1f[feats.mod1.5_l1, !sample])))
table(prediction = pred.mod1.5_l1f, real = vsd.ori$diagnosis2[!sample] )


## TCGA
tab1.4 <- read.table("results/TCGA_gexp/v1.4/model_features/dense.tsv", header = TRUE)
sum(colMeans(tab1.4 == 0) < 1)
mean(colMeans(tab1.4 == 0) < 1)
tab1.4f <- tab1.4[, colMeans(tab1.4 == 0) < 1]

tab1.5 <- read.table("results/TCGA_gexp/v1.5/model_features/dense.tsv", header = TRUE)
sum(colMeans(tab1.5 == 0) < 1)
mean(colMeans(tab1.5 == 0) < 1)
tab1.5f <- tab1.5[, colMeans(tab1.5 == 0) < 1]

tab_auto1.4 <- read.table("results/TCGA_gexp/autoencod_v1.4/model_features/dense.tsv", header = TRUE)
sum(colMeans(tab_auto1.4 == 0) < 1)
mean(colMeans(tab_auto1.4 == 0) < 1)
tab_auto1.4f <- tab_auto1.4[, colMeans(tab_auto1.4 == 0) < 1]

tab_auto1.5 <- read.table("results/TCGA_gexp/autoencod_v1.5/model_features/dense.tsv", header = TRUE)
sum(colMeans(tab_auto1.5 == 0) < 1)
mean(colMeans(tab_auto1.5 == 0) < 1)
tab_auto1.5f <- tab_auto1.5[, colMeans(tab_auto1.5 == 0) < 1]


cor_4.5 <- cor(tab1.4f, tab1.5f)
cor_auto4 <- cor(tab1.4f, tab_auto1.4f )
