#'#################################################################################
#'#################################################################################
#' Explore hsa00430
#'#################################################################################
#'#################################################################################


## Load libraries
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)
library(pheatmap)
library(org.Hs.eg.db)
library(corrplot)
library(qrnn)

genes <- read.table("./results/GTEx_coding/input_genes.txt")$V1
path.map <- read.table("results/GTEx_coding/go_kegg_final_gene_map.tsv", header = TRUE)

gtex.feat <- read.table("results/GTEx_coding/paths_filt2_full_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(gtex.feat) <- paths.vec

weights <- h5read("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/model_weights.h5","weights_paths")
rownames(weights) <- paths.vec
colnames(weights) <- genes

gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")

weights_list <- lapply(letters[1:5], function(submod){
  w <- h5read(paste0("results/GTEx_coding/paths_filt2_full_v3.6", submod, "/model_trained/model_weights.h5"),"weights_paths")
  rownames(w) <- paths.vec
  colnames(w) <- genes
  w
})

all_gtex <- h5read("results/GTEx/all_reshaped_standardized.h5","methy")
rownames(all_gtex) <- genes

## Select hsa00430 weights
path_w <- weights["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol]

path_df <- data.frame(ensembl = names(path_w), weights = path_w)
path_cor <- cor(t(data.matrix(assay(gtex.vst[path_df$ensembl, ]))), gtex.feat[, "path:hsa00430"])

path_df$cor <- path_cor
path_df$Symbol <- mapIds(org.Hs.eg.db, path_df$ensembl , keytype= "ENSEMBL", column="SYMBOL")


gene_cors <- cor(t(data.matrix(assay(gtex.vst[path_df$ensembl, ]))))
colnames(gene_cors) <- rownames(gene_cors) <- path_df$Symbol
corrplot(gene_cors, method = "number", order = "hclust")

pc_path <- prcomp(t(data.matrix(assay(gtex.vst[path_df$ensembl, ]))), scale = TRUE)

plot(pc_path$x[, 1], gtex.feat[, "path:hsa00430"] )


### Get weights of replicates
path_w_mat <- sapply(weights_list, function(m) m["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol])
path_w_mat <- cbind(path_w, path_w_mat )
colnames(path_w_mat) <- c("main", letters[1:5])

ini_w <-  h5read("results/GTEx_coding/paths_filt2_full_v3.6/hsa00430_ini_weights.h5","weights_paths")

path_w_comb <- cbind(path_w_mat, t(ini_w))
colnames(path_w_comb) <- c(colnames(path_w_mat), 1:10)
path_w_comb[, c(3, 5, 8, 12, 14)] <- -path_w_comb[, c(3, 5, 8, 12, 14)]
rownames(path_w_comb) <- path_df$Symbol

heatmap(path_w_comb, scale = "column", col = cm.colors(256))


corrplot(cor(t(path_w_comb)), method = "number", order = "hclust")


path_w_df <- t(path_w_comb)/colSums(abs(path_w_comb))
path_w_df <- data.frame(path_w_df) %>%
  mutate(dataset = colnames(path_w_comb)) %>%
  gather(Gene, value, 1:16) %>%
  mutate(training = ifelse(dataset %in% 1:10, "Pre", "Full"),
        training = factor(training, levels = c("Pre", "Full")))

ggplot(path_w_df, aes(x = Gene, color = training, y = value)) +
  geom_boxplot() +
  theme_bw()

pre_path <- t(all_gtex[path_df$ensembl, ]) %*% path_w_comb[, 8]
plot(elu(pre_path),  gtex.feat[, "path:hsa00430"])
plot(pre_path,  pc_path$x[, 1])

cor(elu(pre_path),  gtex.feat[, "path:hsa00430"])
cor(pre_path,  pc_path$x[, 1])

path_genes <- data.frame(t(assay(gtex.vst[path_df$ensembl, ])))
path_genes$tissue <- gtex.vst$smts

path_df$tissue_var <- sapply(path_df$ensembl, function(gene) {
   summary(lm(formula(paste(gene,  "~ tissue")), path_genes))$adj.r.squared
 })






## Plot pathway values
df_path <- data.frame(path_pre = elu(pre_path), path_full =  gtex.feat[, "path:hsa00430"],
                      pc = pc_path$x[, 1], tissue = gtex.vst$smts)
df_path_filt <- subset(df_path, !tissue %in% c("", "Fallopian Tube", "Bladder", "Cervix Uteri", "Kidney"))
df_path_filt %>%
  ggplot(aes(x = tissue, y = path_pre)) +
  geom_boxplot() +
  theme_bw()



#
df_path_filt %>%
  ggplot(aes(x = tissue, y = path_full)) +
  geom_boxplot() +
  theme_bw()

#
df_path_filt %>%
  ggplot(aes(x = tissue, y = pc)) +
  geom_boxplot() +
  theme_bw()

df_path_filt %>%
  mutate(col = ifelse(!tissue %in% c("Liver", "Brain", "Bone Marrow", "Colon", "Prostate", "Stomach"), "Other tissue", tissue)) %>%
    ggplot(aes(x = path_pre, y = path_full, col = col)) +
    geom_point() +
    theme_bw()

w_mix <- path_w_comb[, 8]
w_mix[c("FMO1", "CDO1", "FMO2", "BAAT")] <- path_w_comb[c("FMO1", "CDO1",  "FMO2",  "BAAT"), 1]
path_mix <- t(all_gtex[path_df$ensembl, ]) %*% w_mix

df_path %>%
  mutate(col = ifelse(!tissue %in% c("Liver", "Brain", "Bone Marrow", "Colon", "Prostate", "Stomach"), "Other tissue", tissue),
          mix = elu(path_mix)) %>%
      ggplot(aes(x = path_pre, y = mix, col = col)) +
      geom_point() +
      theme_bw()


#
df_path %>%
  mutate(col = ifelse(!tissue %in%  c("Liver", "Brain", "Bone Marrow", "Colon", "Prostate", "Stomach"), "Other tissue", tissue),
          mix = elu(path_mix)) %>%
      ggplot(aes(x = mix, y = path_full, col = col)) +
      geom_point() +
      theme_bw()

#
w_mix <- path_w_comb[, 8]
w_mix[c("GAD1", "GAD2", "ADO")] <- path_w_comb[c("GAD1", "GAD2", "ADO"), 1]
path_mix <- t(all_gtex[path_df$ensembl, ]) %*% w_mix

df_path %>%
  mutate(col = ifelse(!tissue %in% c("Liver", "Brain", "Bone Marrow", "Colon", "Prostate", "Stomach"), "Other tissue", tissue),
          mix = elu(path_mix)) %>%
      ggplot(aes(x = path_pre, y = mix, col = col)) +
      geom_point() +
      theme_bw()


#
df_path %>%
  mutate(col = ifelse(!tissue %in% c("Liver", "Brain", "Bone Marrow", "Colon", "Prostate", "Stomach"), "Other tissue", tissue),
          mix = elu(path_mix)) %>%
      ggplot(aes(x = mix, y = path_full, col = col)) +
      geom_point() +
      theme_bw()

png("figures/diff_path_activation_hs00430.png", width = 2000)
df_path_filt %>%
  mutate(path_full  = path_full * 0.930495 +  0.106673 ) %>%
  gather(train, value, 1:2) %>%
  mutate(train = factor(train, levels = c("path_pre", "path_full"))) %>%
  ggplot(aes(x = tissue, y = value, col = train)) +
  geom_boxplot() +
  theme_bw()
dev.off()


## Select hsa00100 weights
path_w2 <- weights["path:hsa00100", subset(path.map, PathwayID == "path:hsa00100")$Symbol]

path_df2 <- data.frame(ensembl = names(path_w2), weights = path_w2)
path_cor2 <- cor(t(data.matrix(assay(gtex.vst[path_df2$ensembl, ]))), gtex.feat[, "path:hsa00100"])

path_df2$cor <- path_cor2
path_df2$Symbol <- mapIds(org.Hs.eg.db, path_df2$ensembl , keytype= "ENSEMBL", column="SYMBOL")


gene_cors2 <- cor(t(data.matrix(assay(gtex.vst[path_df2$ensembl, ]))))
colnames(gene_cors2) <- rownames(gene_cors2) <- path_df2$Symbol
corrplot(gene_cors2, method = "number", order = "hclust")

pc_path2 <- prcomp(t(data.matrix(assay(gtex.vst[path_df2$ensembl, ]))), scale = TRUE)

plot(pc_path2$x[, 1], gtex.feat[, "path:hsa00100"] )



## Select hsa00220 weights
path_w3 <- weights["path:hsa00220", subset(path.map, PathwayID == "path:hsa00220")$Symbol]

path_df3 <- data.frame(ensembl = names(path_w3), weights = path_w3)
path_cor3 <- cor(t(data.matrix(assay(gtex.vst[path_df3$ensembl, ]))), gtex.feat[, "path:hsa00220"])

path_df3$cor <- path_cor3
path_df3$Symbol <- mapIds(org.Hs.eg.db, path_df3$ensembl , keytype= "ENSEMBL", column="SYMBOL")


gene_cors3 <- cor(t(data.matrix(assay(gtex.vst[path_df3$ensembl, ]))))
colnames(gene_cors3) <- rownames(gene_cors3) <- path_df3$Symbol
corrplot(gene_cors3, method = "number", order = "hclust")

pc_path3 <- prcomp(t(data.matrix(assay(gtex.vst[path_df3$ensembl, ]))), scale = TRUE)

plot(pc_path3$x[, 1], gtex.feat[, "path:hsa00220"] )
