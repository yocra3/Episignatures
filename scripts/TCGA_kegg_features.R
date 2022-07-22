#'#################################################################################
#'#################################################################################
#' Explore GSE57945 features from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(HDF5Array)
library(BiocParallel)
library(rjson)
library(cowplot)
library(ggcorrplot)

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

read_training <- function(path, name){

  mat <- read.table(path, row.names = 1) %>% t()
  df <- mat[, c(1, 3)] %>%
      data.frame() %>%
      mutate(epoch = seq_len(nrow(mat))) %>%
      gather(measure, mse, 1:2) %>%
      mutate(model = name)
}

## Compare trainings
all.paths <- read_training("results/GTEx_coding/paths_all_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Gene Sets (Initial)")
base <- read_training("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Gene Set")
base2 <- read_training("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding_training_evaluation.txt", "Whole training")

drop <- read_training("results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7/model_trained/GTEx_coding_training_evaluation.txt", "Step 1 + step 3 + dropout")
drop.full <- read_training("results/GTEx_coding/paths_filt2_full_drop_prime_v3.9/model_trained/GTEx_coding_training_evaluation.txt", "Whole training + dropout")

pre <- read_training("results/GTEx_coding/paths_filt2_full_predense_v6.2/model_trained/GTEx_coding_training_evaluation.txt", "Dense + Gene Set")
post <- read_training("results/GTEx_coding/paths_filt2_full_postdense_v4.3/model_trained/GTEx_coding_training_evaluation.txt", "Gene Set + Dense")
pre_post <- read_training("results/GTEx_coding/paths_filt2_full_prepostdense_v5.3/model_trained/GTEx_coding_training_evaluation.txt", "Dense + Gene Set + Dense")
auto <- read_training("results/GTEx_coding/autoencod_v2.3/model_trained/GTEx_coding_training_evaluation.txt", "Dense autoencoder")


df.models <- Reduce(rbind, list(base, post, pre, pre_post)) %>%
  mutate(model = factor(model, levels = c("Gene Set", "Gene Set + Dense", "Dense + Gene Set", "Dense + Gene Set + Dense" )),
          dataset = ifelse(measure == "loss", "Training", "Validation"))


#
png("figures/TCGA.pathways.trainingeval_models.png", width = 900)
df.models %>%
filter(epoch > 1) %>%
  ggplot(aes(x = epoch, y = mse, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  facet_grid(~ model) +
  theme_bw() +
  scale_color_discrete(name = "") +
  scale_y_continuous(name = "MSE")
dev.off()

df.train <- Reduce(rbind, list(base2, drop, drop.full, auto)) %>%
  mutate(model = factor(model, levels = c("Whole training", "Whole training + dropout", "Step 1 + step 3 + dropout", "Dense autoencoder")),
          dataset = ifelse(measure == "loss", "Training", "Validation"))


png("figures/TCGA.pathways.trainingeval_training.png", width = 700, height = 300)
df.train  %>%
filter(epoch > 1) %>%
ggplot(aes(x = epoch, y = mse, color = dataset, group = dataset)) +
  geom_point() +
  geom_line() +
  facet_grid( ~ model) +
  theme_bw() +
  scale_color_discrete(name = "") +
  scale_y_continuous(name = "MSE")
dev.off()

#
# png("figures/TCGA.kegg.training.epoch10.png", width = 900)
# filter(all.vals, epoch > 10) %>%
# ggplot(aes(x = epoch, y = mse, color = measure, group = measure)) +
#   geom_point() +
#   geom_line() +
#   facet_grid(type ~ model) +
#   theme_bw()
# dev.off()
#
# png("figures/TCGA.kegg.training.epoch10b.png", width = 600)
# filter(all.vals, epoch > 10) %>%
# ggplot(aes(x = epoch, y = mse, color = measure, group = measure)) +
#   geom_point() +
#   geom_line() +
#   facet_grid(model ~ type) +
#   theme_bw()
# dev.off()

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


#

## Compute mse
ori.train <- h5read("results/GTEx/all_reshaped_standardized.h5","methy")
ori.prad <- h5read("results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5","methy")


base <- h5read("results/GTEx_coding/paths_filt2_full_v3.11/model_features/autoencoder_output.h5", "auto")
prad <- h5read("results/GTEx_coding_PRAD/paths_filt2_full_v3.11/model_features/autoencoder_output.h5", "auto")

train_ind <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/test_indices.csv")$V1 + 1


ori.ibd <- h5read("results/SRP042228/assay_reshaped_coding_std_gse.h5", "methy")
ibd <- h5read("results/SRP042228/paths_filt2_full_v3.11/model_features/autoencoder_output.h5", "auto")


# mat_list <- list(base, post, pre, pre_post)
# names(mat_list) <- c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")
vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")
vst_prad <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad")
vst_ibd <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm_TCGA_codingGenes_")

labels <- as.character(read.table("./results/GTEx/individuals_labels.txt", sep = "\t")$V1)
labels.prad <- as.character(read.table("./results/TCGA_gexp_combat_coding/individuals_labels.txt")$V1)

genes <- read.table("./results/GTEx_coding/input_genes.txt")
rownames(ori.ibd) <- rownames(ibd) <- rownames(ori.prad)  <- rownames(base) <- rownames(ori.train) <- rownames(prad) <- as.character(genes$V1)

# phenos <- TCGAbiolinks:::get_IDs(vst)
phenos.prad <- TCGAbiolinks:::get_IDs(vst_prad)

ori.val <- ori.train[, train_ind]
base.val <- base[, train_ind]


paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])

kegg.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
kegg.N <- table(kegg.map$PathwayID)


kegg.df.com <- subset(kegg.df, pathID %in% paths.vec)
kegg.genes.N <- kegg.map %>%
  group_by(Symbol) %>%
  summarize(N = n())


## MSE
# mean.mse.ori <- mean((ori.val - base.val)**2)
#
# # median.mses <- sapply(mat_list, function(x) median((x - ori)**2))
# genes.mse <- rowMeans((ori.val - base.val)**2)
# genes.cors <- sapply(seq_len(nrow(ori.val)), function(i) cor(ori.val[i, ], base.val[i, ], method = "spearman"))
# names(genes.cors) <- rownames(ori.val)
#
genes.cors2 <- sapply(seq_len(nrow(ori.val)), function(i) cor(ori.val[i, ], base.val[i, ], method = "pearson"))
names(genes.cors2) <- rownames(ori.val)

#
# mean.mse.prad <- mean((ori.prad - prad)**2)
# genes.mse.prad <- rowMeans((ori.prad - prad)**2)
# genes.cors.prad <- sapply(seq_len(nrow(ori.prad)), function(i) cor(ori.prad[i, ], prad[i, ], method = "spearman"))
# names(genes.cors.prad) <- rownames(ori.prad)

genes.cors.prad2 <- sapply(seq_len(nrow(ori.prad)), function(i) cor(ori.prad[i, ], prad[i, ], method = "pearson"))
names(genes.cors.prad2) <- rownames(ori.prad)

# genes.mse.ibd <- rowMeans((ori.ibd - ibd)**2)
# genes.cors.ibd <- sapply(seq_len(nrow(ori.ibd)), function(i) cor(ori.ibd[i, ], ibd[i, ], method = "spearman"))
# names(genes.cors.ibd) <- rownames(ori.ibd)
genes.cors.ibd2 <- sapply(seq_len(nrow(ori.ibd)), function(i) cor(ori.ibd[i, ], ibd[i, ], method = "pearson"))
names(genes.cors.ibd2) <- rownames(ori.ibd)


## R2
# r2.ori <- 1 - sum((ori.val - base.val)**2)/sum(( ori.val - rowMeans(ori.val))**2)
# r2.prad <- 1 - sum((ori.prad - prad)**2)/sum((ori.prad - rowMeans(ori.prad))**2)
#
# genes.r2 <- 1 - rowSums((ori.val - base.val)**2)/rowSums(( ori.val - rowMeans(ori.val))**2)
# genes.r2.prad <- 1 - rowSums((ori.prad - prad)**2)/rowSums((ori.prad - rowMeans(ori.prad))**2)
# genes.r2.ibd <- 1 - rowSums((ori.ibd - ibd)**2)/rowSums((ori.ibd - rowMeans(ori.ibd))**2)


# df.cors_mse <- tibble(MSE = c(genes.mse, genes.mse.prad, genes.mse.ibd), correlation = c(genes.cors, genes.cors.prad, genes.cors.ibd),
#   r2 = c(genes.r2, genes.r2.prad, genes.r2.ibd), r2_lineal = c(genes.cors2, genes.cors.prad2, genes.cors.ibd2)**2,
df.cors_mse <- tibble( r2_lineal = c(genes.cors2, genes.cors.prad2, genes.cors.ibd2)**2,
  Symbol = c(names(genes.cors2), names(genes.cors.prad2),names(genes.cors.ibd2)),
  Dataset = rep(c("GTEx validation", "PRAD", "IBD"), c(length(genes.cors2), length(genes.cors.prad2), length(genes.cors.ibd2))))   %>%
  left_join(kegg.genes.N, by = "Symbol") %>%
  mutate(N = ifelse(is.na(N), 0, N),
          path = ifelse(N == 0, "out", "in"),
          Dataset = factor(Dataset, levels = c("GTEx validation", "PRAD", "IBD")))


#
# png("figures/mse_cor_Npathway.png")
# df.cors_mse  %>%
#   mutate(MSE = ifelse(MSE > 2.4, 2.4, MSE),
#         class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
#          Class = factor(class, levels = c("0", "1", "2", "3-5", "6+"))) %>%
#   gather(Measure, Value, 1:2) %>%
#   ggplot(aes(x = Class, y = Value)) +
#   geom_boxplot() +
#   xlab("N pathways per gene") +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scales  = "free_y")
# dev.off()
#
# png("figures/r2_Npathway.png")
# df.cors_mse  %>%
#   mutate(r2 = ifelse(r2 < -1, -1, r2),
#         class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
#          Class = factor(class, levels = c("0", "1", "2", "3-5", "6+"))) %>%
#   gather(Measure, Value, 3:4) %>%
#   ggplot(aes(x = Class, y = Value)) +
#   geom_boxplot() +
#   xlab("N pathways per gene") +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scales  = "free_y")
# dev.off()

r2_gene_plot <- df.cors_mse  %>%
  mutate(class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
         Class = factor(class, levels = c("0", "1", "2", "3-5", "6+"))) %>%
  ggplot(aes(x = Class, y = r2_lineal)) +
  geom_boxplot() +
  xlab("N pathways per gene") +
  ylab(expression(R^2)) +
  theme_bw() +
  facet_grid( ~ Dataset)

png("figures/r2lineal_Npathway.png", height = 300)
r2_gene_plot
dev.off()


# df.cors_mse %>% group_by(Dataset, path) %>%
#   summarize(m30 = mean(correlation > 0.3, na.rm = T), m50 = mean(correlation > 0.5, na.rm = T),
#   m70 = mean(correlation > 0.7, na.rm = T), m90 = mean(correlation > 0.9, na.rm = T))
#
# #
# df.cors_mse %>% group_by(Dataset, path) %>%
#   summarize(r0 = mean(r2 > 0, na.rm = T), r25 = mean(r2 > 0.25, na.rm = T),
#   r50 = mean(r2 > 0.5, na.rm = T), r75 = mean(r2 > 0.75, na.rm = T))

#
df.cors_mse %>% group_by(Dataset, path) %>%
  summarize(r0 = mean(r2_lineal > 0, na.rm = T), r25 = mean(r2_lineal > 0.25, na.rm = T),
  r50 = mean(r2_lineal > 0.5, na.rm = T), r75 = mean(r2_lineal > 0.75, na.rm = T))
df.cors_mse  %>% group_by(Dataset) %>%
  summarize(r0 = mean(r2_lineal > 0, na.rm = T), r25 = mean(r2_lineal > 0.25, na.rm = T),
  r50 = mean(r2_lineal > 0.5, na.rm = T), r75 = mean(r2_lineal > 0.75, na.rm = T),
  m = median(r2_lineal, na.rm = T))

# summary(lm(MSE ~path, df.cors_mse, subset = Dataset != "PRAD"))
# summary(lm(MSE ~path, df.cors_mse, subset = Dataset == "PRAD"))
summary(lm(r2_lineal ~ Dataset, df.cors_mse))


summary(lm(r2_lineal ~N, df.cors_mse, subset = Dataset == "PRAD"))
summary(lm(r2_lineal ~N, df.cors_mse, subset = Dataset == "GTEx validation"))
summary(lm(r2_lineal ~N, df.cors_mse, subset = Dataset == "IBD"))


# png("figures/mse_category.png", width = 1000)
# df.cors_mse  %>%
#   filter(N > 0) %>%
#   select(-N, -path) %>%
#   gather(Category, N, 5:10) %>%
#   filter(N > 0) %>%
#   gather(measure, Value, 3:4) %>%
#   ggplot(aes(x = Category, y = Value)) +
#   geom_boxplot() +
#   facet_grid(Model ~ measure) +
#   theme_bw()
# dev.off()

# png("figures/mse_gene_type.png", width = 1000)
# mse.df %>%
#   mutate(type = ifelse(gene_type %in% c("protein_coding", "processed_pseudogene", "lncRNA", "unprocessed_pseudogene"), gene_type, "Others"),
#           network = ifelse(N > 0, "Network", "Excluded")) %>%
#   gather(measure, val, 2:3) %>%
#   ggplot(aes(x = type, y = val, col = network)) +
#   geom_boxplot() +
#   geom_hline(yintercept = mean.mse, linetype  = 'dashed') +
#   geom_hline(yintercept = median(genes.mse), linetype  = 'dashed', color = 'blue') +
#   theme_bw() +
#   facet_grid(set ~ measure)
# dev.off()


#
# inds.mse.ori <- colMeans((base.val - ori.val)**2)
# inds.mse.prad <- colMeans((ori.prad - prad)**2)
#
# inds.r2 <- 1 - colSums((ori.val - base.val)**2)/colSums(( ori.val- colMeans(ori.val))**2)
# inds.r2.prad <- 1 - colSums((ori.prad - prad)**2)/colSums((ori.prad - colMeans(ori.prad))**2)

inds.cors2 <- sapply(seq_len(ncol(ori.val)), function(i) cor(ori.val[, i], base.val[, i], method = "pearson"))
inds.cors.prad2 <- sapply(seq_len(ncol(ori.prad)), function(i) cor(ori.prad[, i], prad[, i], method = "pearson"))

inds.cors.ibd2 <- sapply(seq_len(ncol(ori.ibd)), function(i) cor(ori.ibd[, i], ibd[, i], method = "pearson"))

# df.inds_mse <- tibble(mse = c(inds.mse.ori, inds.mse.prad),
#                 r2 = c(inds.r2, inds.r2.prad), r2_lin = c(inds.cors2, inds.cors.prad2)**2,
df.inds_mse <- tibble(r2_lin = c(inds.cors2, inds.cors.prad2)**2,
                  Dataset = rep(c("GTEx validation", "PRAD"), c(length(train_ind), ncol(vst_prad))),
                      barcode = c(colnames(vst)[train_ind], colnames(vst_prad)),
                      tumor = c(labels[train_ind], rep("TCGA-PRAD", ncol(vst_prad)))) %>%
  mutate(Dataset = factor(Dataset, levels = c("GTEx validation", "PRAD"))) %>%
  left_join(phenos.prad, by = "barcode") %>%
  mutate(Condition = ifelse(condition == "normal", "Normal", "Cancer"))

#
tum.pheno <- phenos.prad[, c("condition", "barcode")]
colnames(tum.pheno) <- c("condition", "ID")
tum.pheno$condition <- ifelse(tum.pheno$condition == "normal", "Normal", "Cancer" )

ibd.pheno <- colData(vst_ibd)[, c("diagnosis", "run")]
colnames(ibd.pheno) <- c("condition", "ID")


df.inds_cor <- tibble(r2_lin = c(inds.cors2, inds.cors.prad2, inds.cors.ibd2)**2,
                  Dataset = rep(c("GTEx validation", "PRAD", "IBD"), c(length(train_ind), ncol(vst_prad), ncol(vst_ibd))),
                  ID = c(colnames(vst)[train_ind], colnames(vst_prad), colnames(vst_ibd))) %>%
  mutate(Dataset = factor(Dataset, levels = c("GTEx validation", "PRAD", "IBD"))) %>%
  left_join(rbind(tum.pheno, data.frame(ibd.pheno)), by = "ID") %>%
  mutate(Condition = recode(condition, Control = "Normal"),
        Condition = ifelse(is.na(Condition), "Normal", Condition),
        Condition = factor(Condition, levels = c("Normal", "Cancer",  "CD", "UC")))

#
# png("figures/mse_individuals_validation.png", width = 3000)
# df.inds_mse %>%
#   ggplot(aes(x = tumor, y = mse)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
#
#
# png("figures/mse_individuals_windsor_validation.png", width = 3000)
# df.inds_mse %>%
#   mutate(MSE = ifelse(mse > 2, 2, MSE)) %>%
#   ggplot(aes(x = tumor, y = mse)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
#
#
# png("figures/mse_individuals_condition_validation.png", height = 400)
# df.inds_mse %>%
#   mutate(MSE = ifelse(mse > 2, 2, mse)) %>%
#   ggplot(aes(x = Condition, y = MSE)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
#
#
# png("figures/mse_individuals_seqcenter_validation.png")
# df.inds_mse %>%
#   mutate(MSE = ifelse(mse > 2, 2, mse)) %>%
#   ggplot(aes(x = center, y = MSE)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()


# png("figures/r2_individuals_validation.png", width = 3000)
# df.inds_mse %>%
#   filter(Dataset == "Validation") %>%
#   mutate(Tissue = tumor) %>%
#   ggplot(aes(x = Tissue, y = r2_lin)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scale = "free", space = "free_x")
# dev.off()


png("figures/r2lineal_individuals_validation.png", width = 3000)
df.inds_mse %>%
  filter(Dataset == "GTEx validation") %>%
  mutate(Tissue = tumor) %>%
  ggplot(aes(x = Tissue, y = r2_lin)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(R^2)) +
  xlab("Tissue")
dev.off()

#
# r2_inds <- df.inds_mse %>%
#   ggplot(aes(x = tumor, y = r2_lin)) +
#   geom_boxplot() +
#   theme_bw() +
#   ylab(expression(R^2)) +
#   facet_grid( ~ Dataset, scale = "free")


#

r2_inds_plot <- df.inds_cor %>%
  ggplot(aes(x = Condition, y = r2_lin)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(R^2)) +
  facet_grid( ~ Dataset, scale = "free")


png("figures/r2_individuals_condition_validation.png", height = 300)
r2_inds_plot
dev.off()


# png("figures/r2lineal_panel.png")
# plot_grid(r2_gene_plot, r2_inds_plot, ncol = 1, labels = c("A", "B"))
# dev.off()
#
# png("figures/r2_individuals_seqcenter_validation.png")
# df.inds_mse %>%
#   gather(Measure, Value, 2:3) %>%
#   ggplot(aes(x = center, y = Value)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scale = "free", space = "free_x")
# dev.off()
#
# png("figures/r2lineal_individuals_seqcenter_validation.png", height = 300)
# df.inds_mse %>%
#   ggplot(aes(x = center, y = r2_lin)) +
#   geom_boxplot() +
#   theme_bw() +
#   ylab(expression(R^2)) +
#   xlab("Sequencing Center") +
#   facet_grid( ~ Dataset, scale = "free", space = "free_x")
# dev.off()
#

# df.inds_mse %>% group_by(Dataset, Condition) %>%
#   summarize(r0 = mean(r2_lin > 0, na.rm = T), r25 = mean(r2_lin > 0.25, na.rm = T),
#   r50 = mean(r2_lin > 0.5, na.rm = T), r75 = mean(r2_lin > 0.75, na.rm = T))
df.inds_cor  %>% group_by(Dataset) %>%
  summarize(r0 = mean(r2_lin > 0, na.rm = T), r25 = mean(r2_lin > 0.25, na.rm = T),
  r50 = mean(r2_lin > 0.5, na.rm = T), r75 = mean(r2_lin > 0.75, na.rm = T),
  m = median(r2_lin, na.rm = T) )

#
df.inds_cor  %>% group_by(Dataset, Condition) %>%
  summarize(r0 = mean(r2_lin > 0, na.rm = T), r25 = mean(r2_lin > 0.25, na.rm = T),
  r50 = mean(r2_lin > 0.5, na.rm = T), r75 = mean(r2_lin > 0.75, na.rm = T),
  m = median(r2_lin, na.rm = T))

# summary(lm(r2_lin ~Condition, df.inds_mse, subset = Dataset == "Validation"))
# summary(lm(r2_lin ~tumor, df.inds_mse, subset = Dataset == "Validation"))
#
# df.inds_mse <- df.inds_mse %>% mutate(Center = factor(center, levels = c("07", "01", "13", "31")))
# summary(lm(r2_lin ~ Center, df.inds_mse, subset = Dataset != "PRAD"))
#
# summary(lm(r2_lin ~Condition, df.inds_cor, subset = Dataset == "IBD"))

#
# inds.mse.path.ori <- colMeans((base.val -  ori.val)[rownames(base.val) %in% kegg.map$Symbol,]**2)
# inds.mse.path.prad <- colMeans((ori.prad - prad)[rownames(prad) %in% kegg.map$Symbol,]**2)

# inds.path.r2 <- 1 - colSums((ori.val - base.val)[rownames(base.val) %in% kegg.map$Symbol,]**2)/colSums(( ori.val- colMeans(ori.val))[rownames(ori.val) %in% kegg.map$Symbol,]**2)
# inds.r2.path.prad <- 1 - colSums((ori.prad - prad)[rownames(prad) %in% kegg.map$Symbol,]**2)/colSums((ori.prad - colMeans(ori.prad))[rownames(ori.prad) %in% kegg.map$Symbol,]**2)

inds.path.cors2 <- sapply(seq_len(ncol(ori.val)), function(i) cor(ori.val[rownames(ori.val) %in% kegg.map$Symbol, i], base.val[rownames(base.val) %in% kegg.map$Symbol, i], method = "pearson"))
inds.cors.path.prad2 <- sapply(seq_len(ncol(ori.prad)), function(i) cor(ori.prad[rownames(ori.prad) %in% kegg.map$Symbol, i], prad[rownames(prad) %in% kegg.map$Symbol, i], method = "pearson"))
inds.cors.path.ibd2 <- sapply(seq_len(ncol(ori.ibd)), function(i) cor(ori.ibd[rownames(ori.ibd) %in% kegg.map$Symbol, i], ibd[rownames(ibd) %in% kegg.map$Symbol, i], method = "pearson"))
#
# df.inds_mse_path <- tibble(mse = c(inds.mse.path.ori, inds.mse.path.prad),
#                           r2 = c(inds.path.r2, inds.r2.path.prad),
df.inds_mse_path <- tibble(r2_lin = c(inds.path.cors2 , inds.cors.path.prad2)**2,
                    Dataset = rep(c("GTEx validation", "PRAD"), c(length(train_ind), ncol(vst_prad))),
                      barcode = c(colnames(vst)[train_ind], colnames(vst_prad)),
                      tumor = c(labels[train_ind], rep("TCGA-PRAD", ncol(vst_prad)))) %>%
  mutate(Dataset = factor(Dataset, levels = c("GTEx validation", "PRAD"))) %>%
  left_join(phenos.prad, by = "barcode") %>%
  mutate(Condition = ifelse(condition == "normal", "Normal", "Cancer"))

#
df.inds_cors_path <- tibble(r2_lin = c(inds.path.cors2, inds.cors.path.prad2, inds.cors.path.ibd2)**2,
                  Dataset = rep(c("GTEx validation", "PRAD", "IBD"), c(length(train_ind), ncol(vst_prad), ncol(vst_ibd))),
                  ID = c(colnames(vst)[train_ind], colnames(vst_prad), colnames(vst_ibd))) %>%
  mutate(Dataset = factor(Dataset, levels = c("GTEx validation", "PRAD", "IBD"))) %>%
  left_join(rbind(tum.pheno, data.frame(ibd.pheno)), by = "ID") %>%
  mutate(Condition = recode(condition, Control = "Normal"),
         Condition = ifelse(is.na(Condition), "Normal", Condition),
           Condition = factor(Condition, levels = c("Normal", "Cancer",  "CD", "UC")))

#
# png("figures/mse_individuals_path_genes_validation.png", width = 3000)
# df.inds_mse_path %>%
#   ggplot(aes(x = tumor, y = mse)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
# png("figures/mse_individuals_path_condition_validation.png", height = 400)
# df.inds_mse_path %>%
#   mutate(MSE = ifelse(mse > 2, 2, mse)) %>%
#   ggplot(aes(x = Condition, y = MSE)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
#
#
# png("figures/mse_individuals_path_seqcenter_validation.png")
# df.inds_mse_path %>%
#   mutate(MSE = ifelse(mse > 2, 2, mse)) %>%
#   ggplot(aes(x = center, y = MSE)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ Dataset, scale = "free_x", space = "free_x")
# dev.off()
#
#
#
#
#
# png("figures/r2_individuals_path_genes_validation.png", width = 3000)
# df.inds_mse_path %>%
#   gather(Measure, Value, 2:3) %>%
#   ggplot(aes(x = tumor, y = Value)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scale = "free", space = "free_x")
# dev.off()
#
# png("figures/r2_individuals_path_condition_validation.png", height = 400)
# df.inds_mse_path %>%
#   gather(Measure, Value, 2:3) %>%
#   ggplot(aes(x = Condition, y = Value)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scale = "free", space = "free_x")
# dev.off()
#
#
#
# png("figures/r2_individuals_path_seqcenter_validation.png")
# df.inds_mse_path %>%
#   gather(Measure, Value, 2:3) %>%
#   ggplot(aes(x = center, y = Value)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(Measure ~ Dataset, scale = "free", space = "free_x")
# dev.off()

png("figures/r2lineal_individuals_pathgenes_condition.png", height = 300)
 df.inds_cors_path %>%
  ggplot(aes(x = Condition, y = r2_lin)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(R^2)) +
  facet_grid( ~ Dataset, scale = "free")
dev.off()


png("figures/r2lineal_individuals_pathgenes_validation.png", width = 3000)
df.inds_mse_path %>%
  filter(Dataset == "GTEx validation") %>%
  mutate(Tissue = tumor) %>%
  ggplot(aes(x = Tissue, y = r2_lin)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(R^2)) +
  xlab("Tissue")
dev.off()

#
# png("figures/r2lineal_individuals_pathgenes_seqcenter_validation.png", height = 300)
# df.inds_mse_path %>%
#   ggplot(aes(x = center, y = r2_lin)) +
#   geom_boxplot() +
#   theme_bw() +
#   ylab(expression(R^2)) +
#   xlab("Sequencing Center") +
#   facet_grid( ~ Dataset, scale = "free", space = "free_x")
# dev.off()
#

df.inds_cors_path %>% group_by(Dataset, Condition) %>%
  summarize(r0 = mean(r2_lin > 0, na.rm = T), r25 = mean(r2_lin > 0.25, na.rm = T),
  r50 = mean(r2_lin > 0.5, na.rm = T), r75 = mean(r2_lin > 0.75, na.rm = T),  m = median(r2_lin, na.rm = T))

df.inds_cors_path  %>% group_by(Dataset) %>%
  summarize(r0 = mean(r2_lin > 0, na.rm = T), r25 = mean(r2_lin > 0.25, na.rm = T),
  r50 = mean(r2_lin > 0.5, na.rm = T), r75 = mean(r2_lin > 0.75, na.rm = T),  m = median(r2_lin, na.rm = T))


summary(lm(r2_lin ~Condition, df.inds_mse_path, subset = Dataset != "PRAD"))
summary(lm(r2_lin ~tumor, df.inds_mse_path, subset = Dataset != "PRAD"))
summary(lm(r2_lin ~Condition, df.inds_cors_path, subset = Dataset == "IBD"))

df.inds_mse_path <- df.inds_mse_path %>% mutate(Center = factor(center, levels = c("07", "01", "13", "31")))
summary(lm(r2_lin ~ Center, df.inds_mse_path, subset = Dataset != "PRAD"))


## Correlation between PCs
pc.prad <- prcomp(t(prad) )
pc.prad.ori <- prcomp(t(ori.prad))

prad.f <- read.table("results/GTEx_coding_PRAD/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
pc.prad.feat <-  prcomp(prad.f)

df_pcs <- data.frame(PC1 = c(pc.prad.ori$x[, 1], pc.prad.feat$x[, 1], pc.prad$x[, 1]),
                     PC2 = c(pc.prad.ori$x[, 2], -pc.prad.feat$x[, 2], pc.prad$x[, 2]),
                      dataset = rep(c("Original gene expression", "Gene set activities", "Reconstructed gene expression"), c(nrow(pc.prad.ori$x), nrow(pc.prad.feat$x), nrow(pc.prad$x))),
                    Sample = rep(factor(vst_prad$sample_type), 3)) %>%
                    mutate(dataset = factor(dataset, levels = c("Original gene expression", "Gene set activities", "Reconstructed gene expression")),
                          Sample = ifelse(Sample == "Solid Tissue Normal", "Normal", "Cancer"),
                          Sample = factor(Sample, levels = c("Normal", "Cancer")))

png("figures/PRAD_PCS_comparative.png", width = 1200)
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("green", "red")) +
  facet_wrap(~ dataset, scales = "free")
dev.off()

summary(lm(PC1 ~ Sample, df_pcs, subset = dataset == "Original gene expression"))
summary(lm(PC1 ~ Sample, df_pcs, subset = dataset == "Gene set activities"))
summary(lm(PC1 ~ Sample, df_pcs, subset = dataset == "Reconstructed gene expression"))

pcs_model <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("green", "red")) +
ggtitle("PRAD dataset") +
  theme(plot.title = element_text(hjust = 0.5))


pc.ibd <- prcomp(t(ibd) )
pc.ibd.ori <- prcomp(t(ori.ibd))

ibd.f <- read.table("results/SRP042228/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
pc.ibd.feat <-  prcomp(ibd.f)

df_pcs_ibd <- data.frame(PC1 = c(pc.ibd.ori$x[, 1], -pc.ibd.feat$x[, 1], pc.ibd$x[, 1]),
                     PC2 = c(pc.ibd.ori$x[, 2], pc.ibd.feat$x[, 2], pc.ibd$x[, 2]),
                      dataset = rep(c("Original gene expression", "Gene set activities", "Reconstructed gene expression"), c(nrow(pc.ibd.ori$x), nrow(pc.ibd.feat$x), nrow(pc.ibd$x))),
                    Sample = rep(factor(vst_ibd$diagnosis), 3)) %>%
                    mutate(dataset = factor(dataset, levels = c("Original gene expression", "Gene set activities", "Reconstructed gene expression")),
                            Sample = recode(Sample, Control = "Normal"),
                            Sample = factor(Sample, levels = c("Normal", "CD", "UC")))

#
pcs_model_ibd <- ggplot(df_pcs_ibd, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("green", "blue", "orange")) +
  ggtitle("IBD dataset") +
  theme(plot.title = element_text(hjust = 0.5))


png("figures/IBD_PCS_comparative.png", width = 1200)
pcs_model_ibd
dev.off()



vars <- cumsum(pc.prad.ori$sdev**2)/sum(pc.prad.ori$sdev**2)
cors <- diag(cor(pc.prad.ori$x, pc.prad.feat$x))

vars_ibd <- cumsum(pc.ibd.ori$sdev**2)/sum(pc.ibd.ori$sdev**2)
cors_ibd <- diag(cor(pc.ibd.ori$x, pc.ibd.feat$x))


# png("figures/PRAD_cor_pcs_ori_model.png", width = 1000)
# tibble(Correlation = abs(c(cors, cors_ibd)), PC = c(seq_len(length(cors)), seq_len(length(cors_ibd))), Variance = c(vars, vars_ibd),
#       Dataset = rep(c("PRAD", "IBD"), c(length(cors), length(cors_ibd)))) %>%
#   filter(PC < 100) %>%
#   gather(Measure, Value, c(1, 3)) %>%
#   mutate(Measure = recode(Measure, Variance = "Cumulative Variance (%)"),
#         Dataset = factor(Dataset, levels = c("PRAD", "IBD"))) %>%
#   ggplot(aes(x = PC, y = Value, group = Measure, color = Measure)) +
#   geom_line() +
#   facet_wrap(~ Dataset) +
#   theme_bw()
# dev.off()

cors_mat <- cor(pc.prad.ori$x, pc.prad.feat$x)[1:20, 1:20]
colMaxs(abs(cors_mat))

corplot_prad <- ggcorrplot(cors_mat, method = "circle", hc.order = FALSE,
      title = "PRAD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "Gene set activities") +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 20),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 20))

png("figures/PRAD_cor_pcs_ori_model_corrplot.png", width = 700)
corplot_prad
dev.off()

cors_mat_ibd <- cor(pc.ibd.ori$x, pc.ibd.feat$x)[1:20, 1:20]
colMaxs(abs(cors_mat_ibd)[1:10, 1:10])

corplot_ibd <- ggcorrplot(cors_mat_ibd, method = "circle", hc.order = FALSE,
      title = "IBD") +
  scale_x_discrete("Original gene expression") +
  scale_y_discrete(name = "Gene set activities") +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
   axis.title.x = element_text(angle = 0, color = 'grey20', size = 20),
   axis.title.y = element_text(angle = 90, color = 'grey20', size = 20))
png("figures/IBD_cor_pcs_ori_model_corrplot.png", width = 700)
corplot_ibd
dev.off()

png("figures/evaluationTCGAPRAD_panel.png", width = 900, height = 1200)
plot_grid(pcs_model, pcs_model_ibd, plot_grid(corplot_prad, corplot_ibd, nrow = 1, labels = c("C", "D")),  ncol = 1, labels = c("A", "B", ""))
dev.off()
