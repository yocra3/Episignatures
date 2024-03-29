#'#################################################################################
#'#################################################################################
#' Explore hsa00430
#'#################################################################################
#'#################################################################################

# Load libraries
library(tidyverse)
library(cowplot)
library(HDF5Array)
library(SummarizedExperiment)
library(pheatmap)
library(org.Hs.eg.db)
library(ggcorrplot)
library(qrnn)

genes <- read.table("results/GTEx_coding/input_genes.txt")$V1
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)

gtex.feat <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(gtex.feat) <- paths.vec

weights <- h5read("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/model_weights.h5","weights_paths")
rownames(weights) <- paths.vec
colnames(weights) <- genes

weights_pre <- h5read("results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/model_weights.h5","weights_paths")
rownames(weights_pre) <- paths.vec
colnames(weights_pre) <- genes

gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")

weights_list <- lapply(letters[1:5], function(submod){
  w <- h5read(paste0("results/GTEx_coding/paths_filt2_full_v3.11", submod, "/model_trained/model_weights.h5"),"weights_paths")
  rownames(w) <- paths.vec
  colnames(w) <- genes
  w
})

weights_pre_list <- lapply(letters[1:5], function(submod){
  w <- h5read(paste0("results/GTEx_coding/paths_filt2_pre_v3.8", submod, "/model_trained/model_weights.h5"),"weights_paths")
  rownames(w) <- paths.vec
  colnames(w) <- genes
  w
})

all_gtex <- h5read("results/GTEx/all_reshaped_standardized.h5","methy")
rownames(all_gtex) <- genes

## Get matrices of weights
path_w <- weights["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol]
path_w_mat <- sapply(weights_list, function(m) m["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol])
path_w_mat <- cbind(path_w, path_w_mat )
colnames(path_w_mat) <- c("main", letters[1:5])

path_w_pre <- weights_pre["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol]
path_w_pre_mat <- sapply(weights_pre_list, function(m) m["path:hsa00430", subset(path.map, PathwayID == "path:hsa00430")$Symbol])
path_w_pre_mat <- cbind(path_w_pre, path_w_pre_mat )
colnames(path_w_pre_mat) <- 1:6

path_w_comb <- cbind(path_w_mat, path_w_pre_mat)
rownames(path_w_comb)  <- mapIds(org.Hs.eg.db, rownames(path_w_comb)  , keytype= "ENSEMBL", column="SYMBOL")
path_w_comb[, c(2, 4, 5)] <- -path_w_comb[,  c(2, 4, 5)]


png("figures/hsa00430_weights_heatmap.png", res = 300, height = 1500, width = 1500)
myBreaks <- c(seq(-max(abs(path_w_comb)), 0, length.out=ceiling(100/2) + 1),
              seq(max(abs(path_w_comb))/100, max(abs(path_w_comb)), length.out=floor(100/2)))
pheatmap(path_w_comb, breaks = myBreaks)
dev.off()

path_df <- data.frame(ensembl = names(path_w), weights = path_w)
path_df$Symbol <- mapIds(org.Hs.eg.db, path_df$ensembl , keytype= "ENSEMBL", column="SYMBOL")

gene_cors <- cor(t(data.matrix(assay(gtex.vst[path_df$ensembl, ]))))
colnames(gene_cors) <- rownames(gene_cors) <- path_df$Symbol

png("figures/hsa00430_genes_correlation.png", res = 300, height = 1600, width = 1800)
ggcorrplot(gene_cors, method = "circle", hc.order = TRUE)
dev.off()

pc_path <- prcomp(t(data.matrix(assay(gtex.vst[path_df$ensembl, ]))), scale = TRUE)
pre_path <- t(all_gtex[path_df$ensembl, ]) %*% path_w_comb[, "1"]

# pc_plot <- data.frame(PC1 = pc_path$x[, 1], Pathway = pre_path) %>%
#   ggplot(aes(x = Pathway, y = PC1)) + geom_point() +
#   theme_bw() +
#   xlab("Initial gene set activity score") +
#   ylab("PC1 of gene set genes")
# png("figures/hsa00430_step1_pca.png", height = 300, width = 300, res = 300)
# pc_plot
# dev.off()



## Plot pathway values
df_path <- data.frame(path_pre = elu(pre_path), path_full =  gtex.feat[, "path:hsa00430"],
                      pc = pc_path$x[, 1], tissue = gtex.vst$smts)

plot_train <- df_path %>%
  mutate(Tissue = ifelse(!tissue %in% c("Brain"), "Other tissue", tissue)) %>%
    ggplot(aes(x = path_pre, y = path_full, col = Tissue)) +
    geom_point() +
    theme_bw() +
    xlab("Initial gene set activity score") +
    ylab("Final gene set activity score") +
    scale_color_manual(values = c("#fcfc7e", "grey")) +
    geom_abline(slope = 1, linetype = "dashed")
    # scale_y_continuous(limits = c(-1, 3)) +
    # scale_x_continuous(limits = c(-1, 3))

png("figures/hsa00430_step1_step3.png", height = 900, width = 1500, res = 300)
plot_train
dev.off()

# png("figures/hsa00430_scores_panel.png", width = 600, height = 300)
# plot_grid(pc_plot, plot_train, nrow = 1, labels = c("A", "B"), rel_widths = c(2, 3))
# dev.off()


a <- path_w_comb[, c("main", "1")] %>%
  data.frame() %>%
  mutate(Gene = rownames(path_w_comb),
        diff = abs(main - X1)) %>%
  arrange(abs(diff))


plot_weight <- path_w_comb[, c("main", "1")] %>%
  data.frame() %>%
  mutate(Gene = rownames(path_w_comb)) %>%
  filter(Gene %in% c("GAD1", "GAD2", "FMO1")) %>%
  gather(Model, Weight, 1:2) %>%
  mutate(Model = recode(Model, main = "Final", X1 = "Initial"),
          Model = factor(Model, levels = c("Initial", "Final"))) %>%
    ggplot(aes(x = Model, y = abs(Weight), group = Gene, col = Gene)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(name = "Relevance") +
    theme_bw()

png("figures/hsa00430_weights.png", height = 600, width = 900, res = 300)
plot_weight
dev.off()

#
brain_path <- subset(path_df,Symbol  %in% c("GAD1", "GAD2"))
brain_mat <- t(data.matrix(assay(gtex.vst[brain_path$ensembl, ] )))
colnames(brain_mat) <- brain_path$Symbol

plot_genes <- data.frame(brain_mat, Tissue = gtex.vst$smts) %>%
  filter(!Tissue %in% c("", "Fallopian Tube", "Bladder", "Cervix Uteri", "Kidney")) %>%
  mutate(Group = ifelse(!Tissue %in% c("Brain"), "Other tissues", Tissue)) %>%
  gather(Gene, Expression, 1:2) %>%
  ggplot(aes(x = Tissue, y = Expression, fill = Group)) +
  geom_violin() +
  scale_fill_manual(values = c("#f2f2057c", "#bebebe82")) +
  facet_grid(Gene ~ ., scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
png("figures/hsa00430_sel_gene_exprs.png", width = 6000, height = 2000, res = 300)
plot_genes
dev.off()
