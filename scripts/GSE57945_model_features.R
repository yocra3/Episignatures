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
library(sva)


makePCdf <- function(seobj){
  pc <- prcomp(t(assay(seobj)))
  pcdf <- data.frame(pc$x[, 1:10])
  pcdf <- cbind(pcdf, colData(seobj)[, c("diagnosis2", "age", "sex", "age_group")])
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
           show_rownames = TRUE)

}
trainSVM <-  function(seobj){
  mat <- data.matrix(assay(seobj))
  df <- data.frame(pathClass = factor(seobj$diagnosis2), t(mat))
  model_svm <- svm(pathClass ~ ., df)
}

getPathwayCor <- function(path, se_pc, se_path){
  path2 <- gsub(".", ":", path, fixed = TRUE)
  genes <- subset(path_map, PathwayID == path2)$Symbol
  pc <- prcomp(t(assay(se_pc[rownames(se_pc) %in% genes,])), rank. = 1)$x[, 1]
  cor(pc, t(assay(se_path[path,])))

}



## Compute mse
paths <- read.table("results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])

kegg.map <- read.table("results/preprocess/go_kegg_final_gene_map.tsv", header = TRUE)
kegg.N <- table(kegg.map$PathwayID)
kegg.genes.N <- kegg.map %>%
  group_by(Symbol) %>%
  summarize(N = n())

ori <- h5read("results/SRP042228/assay_reshaped_coding_std_gse.h5","methy")


base <- h5read("results/SRP042228/comb_paths3_v3.6/model_features/autoencoder_output.h5", "auto")
auto <- h5read("results/SRP042228/autoencod_v2.3/model_features/autoencoder_output.h5", "auto")

vst <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm")


genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
rownames(ori) <- rownames(base) <- rownames(auto) <- as.character(genes$V1)



## MSE
mean.mse.base <- mean((ori - base)**2)
mean.mse.auto <- mean((ori - auto)**2)

genes.mse.base <- rowMeans((ori - base)**2)
genes.cors.base <- sapply(seq_len(nrow(ori)), function(i) cor(ori[i, ], base[i, ], method = "spearman"))
names(genes.cors.base) <- rownames(ori)

genes.cors.base2 <- sapply(seq_len(nrow(ori)), function(i) cor(ori[i, ], base[i, ], method = "pearson"))

mean.mse.auto <- mean((ori - auto)**2)
genes.mse.auto <- rowMeans((ori - auto)**2)
genes.cors.auto <- sapply(seq_len(nrow(ori)), function(i) cor(ori[i, ], auto[i, ], method = "spearman"))
names(genes.cors.auto) <- rownames(ori)

genes.cors.auto2 <- sapply(seq_len(nrow(ori)), function(i) cor(ori[i, ], auto[i, ], method = "pearson"))

## R2
r2.base <- 1 - sum((ori - base)**2)/sum(( ori - rowMeans(ori))**2)
r2.auto <- 1 - sum((ori - auto)**2)/sum((ori - rowMeans(ori))**2)

genes.r2.base <- 1 - rowSums((ori - base)**2)/rowSums(( ori - rowMeans(ori))**2)
genes.r2.auto <- 1 - rowSums((ori - auto)**2)/rowSums((ori - rowMeans(ori))**2)


df.cors_mse <- tibble(MSE = c(genes.mse.base, genes.mse.auto), correlation = c(genes.cors.auto, genes.cors.auto),
  r2 = c(genes.r2.base, genes.r2.auto), r2_lineal = c(genes.cors.base2, genes.cors.auto2)**2,
  Symbol = c(names(genes.cors.base), names(genes.cors.auto)),
  Dataset = rep(c("Model", "Autoencoder"), c(length(genes.cors.base), length(genes.cors.auto))))   %>%
  left_join(kegg.genes.N, by = "Symbol") %>%
  mutate(N = ifelse(is.na(N), 0, N),
          path = ifelse(N == 0, "out", "in"),
          Dataset = factor(Dataset, levels = c("Model", "Autoencoder")))



#
df.cors_mse  %>%
  mutate(MSE = ifelse(MSE > 2.4, 2.4, MSE),
        class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
         Class = factor(class, levels = c("0", "1", "2", "3-5", "6+"))) %>%
  gather(Measure, Value, 1:2) %>%
  ggplot(aes(x = Class, y = Value)) +
  geom_boxplot() +
  xlab("N pathways per gene") +
  theme_bw() +
  facet_grid(Measure ~ Dataset, scales  = "free_y")


## Load vsd data ####
vsd.ori <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm_TCGA_codingGenes_")
vsd.ori$diagnosis2 <- relevel(factor(vsd.ori$diagnosis2), ref = "Control")
## Create train and test indexes
sample <- sample.split(vsd.ori$diagnosis2, SplitRatio = .80)


## Get features associated with CD
mod <- model.matrix(~diagnosis2 + age + sex, colData(vsd.ori))

lmori <- lmFit(assay(vsd.ori), mod) %>% eBayes()
tab.ori <- topTable(lmori, coef = 2:4, n = Inf)
featsori <- rownames(subset(tab.ori, adj.P.Val  < 0.05))
tab.ori.icd <- topTable(lmori, coef = 3, n = Inf)
tab.ori.ccd <- topTable(lmori, coef = 2, n = Inf)
tab.ori.uc <- topTable(lmori, coef = 4, n = Inf)
# #
# # ## See age
# # tab.ori.age <- topTable(lmori, coef = 5, n = Inf)
# # featsori.age <- rownames(subset(tab.ori.age, adj.P.Val  < 0.05))
# #
#
# ## See sex
# tab.ori.sex <- topTable(lmori, coef = 6, n = Inf)
# featsori.sex <- rownames(subset(tab.ori.sex, adj.P.Val  < 0.05))

## Make PCs
pc.ori <- makePCdf(vsd.ori)
# pc.ori.feats <- makePCdf(vsd.ori[featsori, ])


# sample2 <- sample[pc.ori$PC2 > -50]
# #
# #
# # png("figures/SRP042228.raw.pca.png")
# # makePCplot(pc.ori, "diagnosis2")
# # dev.off()
# #
# # makePCplot(pc.ori, "age_group")
# # makePCplot(pc.ori, "sex")
# #
# # makePCplot(pc.ori.feats, "diagnosis2")
# # makePCplot(pc.ori.feats, "age")
# # makePCplot(pc.ori.feats, "sex")
#
# ## Heatmaps
# #makeHeatmap(vsd.ori)
# makeHeatmap(vsd.ori[featsori, ])

## SVM
svm.ori <- trainSVM(vsd.ori[featsori, sample])
pred.ori <- predict(svm.ori, t(assay(vsd.ori[featsori, !sample])))
table(prediction = pred.ori, real = vsd.ori$diagnosis2[!sample] )

svm.ori.paper <- trainSVM(vsd.ori[featsori, sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])
pred.ori.paper <- predict(svm.ori.paper, t(assay(vsd.ori[featsori, !sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])))
table(prediction = pred.ori.paper, real = vsd.ori$diagnosis2[!sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")] )


#
# vsd.ori2 <- vsd.ori[, pc.ori$PC2 > -50]
# svm.ori2 <- trainSVM(vsd.ori2[featsori, sample2])
# pred.ori2 <- predict(svm.ori2, t(assay(vsd.ori2[featsori, !sample2])))
# table(prediction = pred.ori2, real = vsd.ori2$diagnosis2[!sample2] )
#
# ## UC
# vsd.uc <- vsd.ori[, vsd.ori$diagnosis != "CD"]
# vsd.uc$diagnosis2 <- droplevels(vsd.uc$diagnosis2)
# mod.uc <- model.matrix(~diagnosis2 + age + sex, colData(vsd.uc ))
# lm.uc <- lmFit(assay(vsd.uc), mod.uc) %>% eBayes()
# tab.uc <- topTable(lm.uc, coef = 2, n = Inf)
# feats.uc  <- rownames(subset(tab.uc  , adj.P.Val  < 0.1))
#

### Compute SVA
mod0 <- model.matrix(~ diagnosis2, colData(vsd.ori))
svobj <- sva(data.matrix(assay(vsd.ori)), mod, mod0)

ori.sv <- residuals(lmFit(assay(vsd.ori), svobj$sv), assay(vsd.ori))
vsd.sv <- vsd.ori
assay(vsd.sv) <- ori.sv
pc.sv <- makePCdf(vsd.sv)
cot <- prcomp(t(assay(vsd.sv))()

makePCplot(pc.sv, "diagnosis2")
makePCplot(pc.sv, "age_group")
makePCplot(pc.sv, "sex")

lmori.sv <- lmFit(assay(vsd.ori), cbind(mod, svobj$sv)) %>% eBayes()
tab.ori.sv <- topTable(lmori.sv, coef = 2:4, n = Inf)
featsori.sv <- rownames(subset(tab.ori.sv, adj.P.Val  < 0.05))
tab.sv.ccd <- topTable(lmori.sv, coef = 2, n = Inf)
tab.sv.icd <- topTable(lmori.sv, coef = 3, n = Inf)
write.table(file = "icd_genes.txt", rownames(subset(tab.sv.icd, adj.P.Val < 0.05)), quote = FALSE, row.names = FALSE, col.names = FALSE)

## SVM
svm.ori.sv <- trainSVM(vsd.sv[featsori.sv, sample ])
pred.ori.sv <- predict(svm.ori.sv, t(assay(vsd.sv[featsori.sv, !sample])))
table(prediction = pred.ori.sv, real = vsd.sv$diagnosis2[!sample] )

svm.ori.sv.paper <- trainSVM(vsd.sv[featsori.sv, sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])
pred.ori.sv.paper <- predict(svm.ori.sv.paper, t(assay(vsd.sv[featsori.sv, !sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])))
table(prediction = pred.ori.sv.paper, real = vsd.sv$diagnosis2[!sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")] )


## DNN kegg - pre ####
se.kegg <- readFeatures("results/SRP042228/kegg_filt2_v6.2/model_features/prune_low_magnitude_dense_1.tsv", vsd.ori)
paths <- read.table("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2/model_trained/pathways_names.txt", header = TRUE)
rownames(se.kegg)  <- gsub(":", ".",as.character( paths[, 1]))

## SVM
svm.kegg <- trainSVM(se.kegg[, sample])
pred.kegg <- predict(svm.kegg, t(assay(se.kegg[, !sample])))
table(prediction = pred.kegg, real = vsd.ori$diagnosis2[!sample] )

svm.kegg.paper <- trainSVM(se.kegg[, sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])
pred.kegg.paper <- predict(svm.kegg.paper, t(assay(se.kegg[, !sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])))
table(prediction = pred.kegg.paper , real = vsd.ori$diagnosis2[!sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")] )

svobj.kegg <- sva(data.matrix(assay(se.kegg)), mod, mod0)
kegg.sv <- se.kegg
assay(kegg.sv) <-  residuals(lmFit(assay(se.kegg), svobj.kegg$sv), assay(se.kegg))

svm.kegg.sv <- trainSVM(kegg.sv [, sample])
pred.kegg.sv <- predict(svm.kegg.sv, t(assay(kegg.sv [, !sample])))
table(prediction = pred.kegg.sv, real = vsd.ori$diagnosis2[!sample] )

svm.kegg.sv.paper <- trainSVM(kegg.sv[, sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])
pred.kegg.sv.paper <- predict(svm.kegg.sv.paper, t(assay(kegg.sv[, !sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")])))
table(prediction = pred.kegg.sv.paper, real = kegg.sv$diagnosis2[!sample & vsd.ori$diagnosis2 %in% c("cCD", "UC")] )
#
#
# se.kegg2 <- se.kegg[, pc.ori$PC2 > -50]
# svm.kegg.f <- trainSVM(se.kegg2[, sample2])
# pred.kegg.f <- predict(svm.kegg.f, t(assay(se.kegg2[, !sample2])))
# table(prediction = pred.kegg.f, real = se.kegg2$diagnosis2[!sample2] )

# lm.kegg <- lmFit(assay(se.kegg), mod) %>% eBayes()
# tab.kegg <- topTable(lm.kegg, coef = 2:4, n = Inf)
# feats.kegg  <- rownames(subset(tab.kegg , adj.P.Val  < 0.05))
#

lm.kegg.sv <- lmFit(assay(se.kegg), cbind(mod, svobj$sv)) %>% eBayes()
tab.kegg.sv <- topTable(lm.kegg.sv, coef = 2:4, n = Inf)
feats.kegg.sv  <- rownames(subset(tab.kegg.sv , adj.P.Val  < 0.05))
tab.kegg.icd.sv <- topTable(lm.kegg.sv, coef = 3, n = Inf)


kegg.map <- read.table("results/preprocess/kegg_filt_gene_map.tsv", header = TRUE)
kegg.map$pathID <- gsub(":", ".", kegg.map$PathwayID )

makeHeatmap(vsd.ori[rownames(vsd.ori) %in% subset(kegg.map, pathID == "path.hsa00603")$Symbol,
          vsd.ori$diagnosis2 %in% c("Control", "iCD")])

p <- sapply(rownames(tab.kegg.icd.sv) , function(i) mean(subset(tab.sv.icd, rownames(tab.sv.icd) %in% subset(kegg.map, pathID == i)$Symbol)$adj.P.Val < 0.05))
plot(-log10(tab.kegg.icd.sv$P.Value), p)

tab.kegg.icd.sv$prop_genes <- p
tab.kegg.icd.sv$p_log <- -log10(tab.kegg.icd.sv$P.Value)

summary(robustbase::lmrob(prop_genes ~ p_log, tab.kegg.icd.sv, subset = !is.na(prop_genes) & P.Value < 0.05, k.max = 1000))
summary(lm(prop_genes ~ p_log, tab.kegg.icd.sv, subset = !is.na(prop_genes) & P.Value < 0.5))

png("figures/SRP042228.corr_pathwaysiCD_propGenes.png")
tab.kegg.icd.sv %>%
  filter(P.Value < 0.05) %>%
  ggplot(aes(x = p_log, y = prop_genes)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("-log10 p-value pathway") +
  ylab("Proportion of genes significant")
dev.off()












pc.kegg.sv <- makePCdf(kegg.sv)

makeHeatmap(kegg.sv[, kegg.sv$diagnosis2 != "iCD"])

makePCplot(pc.kegg.sv, "age_group")
makePCplot(pc.kegg.sv, "diagnosis2")

tab.kegg.ccd <- topTable(lm.kegg, coef = 2, n = Inf)
tab.kegg.icd <- topTable(lm.kegg, coef = 3, n = Inf)
tab.kegg.uc <- topTable(lm.kegg, coef = 4, n = Inf)

tab.kegg.icd.sv <- topTable(lm.kegg.sv, coef = 3, n = Inf)

svm.kegg.sv <- trainSVM(kegg.sv[, sample])
pred.kegg.sv <- predict(svm.kegg.sv, t(assay(kegg.sv[, !sample])))
table(prediction = pred.kegg.sv, real = kegg.sv$diagnosis2[!sample] )
table(prediction = predict(svm.kegg.sv, t(assay(kegg.sv[, sample]))), real = kegg.sv$diagnosis2[sample] )

makeHeatmap(se.kegg)

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

feats.kegg.uc  <- rownames(tab.kegg.uc)[1:10]

png("figures/SRP042228.kegg_uc.heatmap.png")
makeHeatmap(se.kegg.uc[feats.kegg.uc , ])
dev.off()

png("figures/SRP042228.kegg_uc.heatmap_all.png")
makeHeatmap(se.kegg.uc)
dev.off()

## Correlations KEGG vs PC1 ####
path_map <- read.table("results/preprocess/kegg_gene_map.tsv", header = TRUE )

### All
path.cores.all <- sapply(rownames(se.kegg), getPathwayCor, se_pc = vsd.ori, se_path = se.kegg)
names(path.cores.all) <- rownames(se.kegg)

png("figures/SRP042228.vsd_genes.kegg.cor.png")
hist(path.cores.all, breaks = 40)
dev.off()


### Control vs UC
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

path.cores.uc <- sapply(rownames(se.kegg.uc), getPathwayCor, vsd.uc, se.kegg.uc)
names(path.cores.uc) <- rownames(se.kegg.uc)
path_n <- sapply(rownames(se.kegg.uc), function(x) sum(gsub(":", ".", path_map$PathwayID) == x))
cor(path_n, abs(path.cores.uc ))

png("figures/SRP042228.uc.vsd_genes.kegg.cor.png")
hist(path.cores.uc, breaks = 40)
dev.off()

path_n2 <- path_n
path_n2[path_n2 > 500] <- 500
png("figures/SRP042228.uc.vsd_genes.kegg.cor_vs_n.png")
plot(path.cores.uc, path_n2)
dev.off()


## DNN hipathia ####
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

## Hipathia ####
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


## Compute correlations Hipathia vs DNN Hipathi #####
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
