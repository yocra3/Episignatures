#'#################################################################################
#'#################################################################################
#' Compare GTEx model with TCGA model
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
library(rjson)
# library(GOfuncR)
library(ggcorrplot)
library(ggrepel)

## Load data
genes <- read.table("./results/GTEx_coding/input_genes.txt")$V1
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)

paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])

weights_gtex <- h5read("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/model_weights.h5","weights_paths")
rownames(weights_gtex) <- paths.vec
colnames(weights_gtex) <- genes

weights_gtex_pre <- h5read("results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/model_weights.h5","weights_paths")
rownames(weights_gtex_pre) <- paths.vec
colnames(weights_gtex_pre) <- genes


weights_tcga <- h5read("results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/model_weights.h5","weights_paths")
rownames(weights_tcga) <- paths.vec
colnames(weights_tcga) <- genes

weights_tcga_pre <- h5read("results/TCGA_coding_all/paths_filt2_pre_v3.8/model_trained/model_weights.h5","weights_paths")
rownames(weights_tcga_pre) <- paths.vec
colnames(weights_tcga_pre) <- genes


gtex.feat <- read.table("results/TCGA_all/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
tcga.feat <- read.table("results/TCGA_all_TCGA/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(gtex.feat) <- colnames(tcga.feat)  <- paths.vec

gtex.feat.pre <- read.table("results/TCGA_all/paths_filt2_pre_v3.8/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
tcga.feat.pre <- read.table("results/TCGA_all_TCGA/paths_filt2_pre_v3.8/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(gtex.feat.pre) <- colnames(tcga.feat.pre) <- paths.vec


gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")
tcga.vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

path.N <- table(path.map$PathwayID) %>%
  data.frame()
#
# ## Load functions
# readPathways <- function(model, sufix, path_name){
#   lapply(sufix, function(i){
#     path <- paste0("results/TCGA_all_TCGA/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
#     tab <- read.table(path, header = TRUE)
#     tab <- data.matrix(tab)
#     colnames(tab) <- path_name
#     tab
#   })
# }
#
# readPathways2 <- function(model, sufix, path_name){
#   lapply(sufix, function(i){
#     path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
#     tab <- read.table(path, header = TRUE)
#     tab <- data.matrix(tab)
#     colnames(tab) <- path_name
#     tab
#   })
# }
#
# pathwayCorr <- function(path_list, col){
#   path_mat <- sapply(path_list, function(x) x[, col])
#   cors <- cor(path_mat)
#   cors[upper.tri(cors)]
# }
#
#
# makeDFsum <- function(cors, mod_name){
#   df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
#     mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
#             model = mod_name) %>%
#             as_tibble()
#
# }
#
# path_vals_tcga <-  readPathways("paths_filt2_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
#
# ## Define replicability
# tcga.cors <- sapply(paths.vec, pathwayCorr, path_list = path_vals_tcga)
# colnames(tcga.cors) <- paths.vec
# df.tcga <- makeDFsum(tcga.cors,"full") %>%
#   left_join(data.frame(path.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
#   mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))
#
# #
# path_vals_gtex <-  readPathways2("paths_filt2_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
#
# ## Define replicability
# gtex.cors <- sapply(paths.vec, pathwayCorr, path_list = path_vals_gtex)
# colnames(gtex.cors) <- paths.vec
# df.gtex <- makeDFsum(gtex.cors,"full") %>%
#   left_join(data.frame(path.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
#   mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

## Compute GTEx - TCGA correlations
weight_cors <- sapply(seq_len(nrow(weights_tcga)), function(i){

  # Gtex
  p_g <- weights_gtex[i, ]
  p_g <- p_g[p_g != 0]

  # TCGA
  p_t <- weights_tcga[i, ]
  p_t <- p_t[p_t != 0]

  cor(p_g, p_t)
})

path_cors <- sapply(seq_len(ncol(gtex.feat)), function(i){

  cor(gtex.feat[, i], tcga.feat[, i])
})

path_cors_pre <- sapply(seq_len(ncol(gtex.feat.pre)), function(i){

  cor(gtex.feat.pre[, i], tcga.feat.pre[, i])
})



path_df <- tibble(path = paths.vec, cor_weight = weight_cors, cor_path = path_cors,
                  cor_path_pre = path_cors_pre)
save(path_df, file = "results/GTEx_coding/paths_filt2_full_v3.11/gtex_tcga_comparative.Rdata")

# ##
# path_df <- mutate(path_df, Category = ifelse(abs(path_cors) > 0.7 & minCor > 0.7 & abs(cor_weight) > 0.7, "High",
#         ifelse(abs(path_cors) < 0.3 & minCor > 0.7,
#           ifelse(abs(cor_path_pre) < 0.3, "Low pre and full",
#               ifelse(abs(cor_path_pre) > 0.7, "Low full", "Other")),
#             "Other")))
#


hist_data_cor <- path_df %>%
  ggplot(aes(x = abs(cor_path))) +
  geom_histogram() +
  ylab("N gene sets") +
  geom_vline(xintercept  = c(0.3, 0.7)) +
  xlab("Correlation between GTEx and TCGA") +
  theme_bw()

png("figures/data_models_cor_hist.png", height = 300)
hist_data_cor
dev.off()


path_df %>%
  summarize(m70 = mean(abs(cor_path) > 0.7), m30 = mean(abs(cor_path) > 0.3))


#
png("figures/data_models_cor_step.png", height = 400)
path_df %>%
  mutate(High = ifelse(path %in% c("GO:0070076", "GO:0010888", "GO:0010528", "GO:0032274"), path, "Other")) %>%
  ggplot(aes(x = abs(cor_path_pre), y = abs(cor_path), color = High)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(name = "Gene sets", values = c("darkgreen", "blue", "orange", "brown", "grey")) +
  xlab("Correlation after step 1") +
  ylab("Correlation after full training")
dev.off()

cor(x = abs(path_df$cor_path_pre), y = abs(path_df$cor_path))

## High pre - high full gene set - GO:0070076
df_h <- data.frame(gtex = weights_gtex_pre["GO:0070076", ], tcga = weights_tcga_pre["GO:0070076", ], gene = colnames(weights_tcga_pre)) %>%
  filter(gtex != 0 & tcga != 0) %>%
  mutate(gene = mapIds(org.Hs.eg.db, gene, keytype= "ENSEMBL", column="SYMBOL"))

plot_hh_pre <- df_h %>%  ggplot(aes(x = gtex, y = tcga)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_h$gtex, df_h$tcga)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0070076 - step 1") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

genes.top <- subset(path.map, PathwayID == "GO:0070076")$Symbol
genes.top_Symbol <- mapIds(org.Hs.eg.db, genes.top, keytype= "ENSEMBL", column="SYMBOL")

gene_cors.top_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.top, ]))))
rownames(gene_cors.top_gtex ) <- colnames(gene_cors.top_gtex ) <- genes.top_Symbol

plot_hh_cor_gtex <-  ggcorrplot(gene_cors.top_gtex, method = "circle", hc.order = TRUE) +
  ggtitle("GO:0070076 - GTEx") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

gene_cors.top_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.top, ]))))
rownames(gene_cors.top_tcga ) <- colnames(gene_cors.top_tcga ) <- genes.top_Symbol

plot_hh_cor_tcga <- ggcorrplot(gene_cors.top_tcga,  method = "circle", hc.order = TRUE) +
  ggtitle("GO:0070076 - TCGA") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))


df_h2 <- data.frame(gtex = -weights_gtex["GO:0070076", ], tcga = weights_tcga["GO:0070076", ], gene = colnames(weights_tcga_pre)) %>%
  filter(gtex != 0 & tcga != 0) %>%
  mutate(gene = mapIds(org.Hs.eg.db, gene, keytype= "ENSEMBL", column="SYMBOL"))
plot_hh_post <-  df_h2 %>%  ggplot(aes(x = gtex, y = tcga)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_h2$gtex, df_h2$tcga)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0070076 - full training") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

png("figures/data_models_high_geneset.png", width = 1100, height = 1100)
plot_grid(plot_hh_pre, plot_hh_post, plot_hh_cor_gtex, plot_hh_cor_tcga, labels = LETTERS[1:4], ncol = 2 )
dev.off()


## Low pre - Low full gene set - GO:0010888
df_path_low2 <- data.frame(gtex_pre = -weights_gtex_pre["GO:0010888", ], tcga_pre = weights_tcga_pre["GO:0010888", ], gene = colnames(weights_tcga_pre),
                          gtex = weights_gtex["GO:0010888", ], tcga = -weights_tcga["GO:0010888", ]) %>%
              filter(gtex != 0 & tcga != 0) %>%
              mutate(gene = mapIds(org.Hs.eg.db, gene, keytype= "ENSEMBL", column="SYMBOL"))

plot_ll_pre <- df_path_low2 %>%  ggplot(aes(x = gtex_pre, y = tcga_pre)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_low2$gtex_pre, df_path_low2$tcga_pre)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0010888 - step 1") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

genes.ll <- subset(path.map, PathwayID == "GO:0010888")$Symbol
genes.ll_Symbol <- mapIds(org.Hs.eg.db, genes.ll, keytype= "ENSEMBL", column="SYMBOL")

gene_cors.ll_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.ll, ]))))
rownames(gene_cors.ll_gtex ) <- colnames(gene_cors.ll_gtex ) <- genes.ll_Symbol

plot_ll_cor_gtex <-  ggcorrplot(gene_cors.ll_gtex, method = "circle", hc.order = TRUE) +
  ggtitle("GO:0010888 - GTEx") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

gene_cors.ll_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.ll, ]))))
rownames(gene_cors.ll_tcga ) <- colnames(gene_cors.ll_tcga ) <- genes.ll_Symbol

plot_ll_cor_tcga <- ggcorrplot(gene_cors.ll_tcga,  method = "circle", hc.order = TRUE) +
  ggtitle("GO:0010888 - TCGA") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_ll_post <-  df_path_low2 %>%  ggplot(aes(x = gtex, y = tcga)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_low2$gtex, df_path_low2$tcga)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0010888 - full training") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

png("figures/data_models_low_geneset.png", width = 1100, height = 1100)
plot_grid(plot_ll_pre, plot_ll_post, plot_ll_cor_gtex, plot_ll_cor_tcga, labels = LETTERS[1:4], ncol = 2 )
dev.off()



## Low pre - high full gene set - GO:0010528
df_path_hilow <- data.frame(gtex_pre = -weights_gtex_pre["GO:0010528", ], tcga_pre = weights_tcga_pre["GO:0010528", ], gene = colnames(weights_tcga_pre),
                          gtex = weights_gtex["GO:0010528", ], tcga = -weights_tcga["GO:0010528", ]) %>%
              filter(gtex != 0 & tcga != 0) %>%
              mutate(gene = mapIds(org.Hs.eg.db, gene, keytype= "ENSEMBL", column="SYMBOL"))

plot_hl_pre <- df_path_hilow %>%  ggplot(aes(x = gtex_pre, y = tcga_pre)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_hilow$gtex_pre, df_path_hilow$tcga_pre)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0010528 - step 1") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

genes.hl <- subset(path.map, PathwayID == "GO:0010528")$Symbol
genes.hl_Symbol <- mapIds(org.Hs.eg.db, genes.hl, keytype= "ENSEMBL", column="SYMBOL")

gene_cors.hl_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.hl, ]))))
rownames(gene_cors.hl_gtex ) <- colnames(gene_cors.hl_gtex ) <- genes.hl_Symbol

plot_hl_cor_gtex <-  ggcorrplot(gene_cors.hl_gtex, method = "circle", hc.order = TRUE) +
  ggtitle("GO:0010528 - GTEx") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

gene_cors.hl_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.hl, ]))))
rownames(gene_cors.hl_tcga ) <- colnames(gene_cors.hl_tcga ) <- genes.hl_Symbol

plot_hl_cor_tcga <- ggcorrplot(gene_cors.hl_tcga,  method = "circle", hc.order = TRUE) +
  ggtitle("GO:0010528 - TCGA") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_hl_post <-  df_path_hilow %>%  ggplot(aes(x = gtex, y = tcga)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_hilow$gtex, df_path_hilow$tcga)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0010528 - full training") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

png("figures/data_models_highlow_geneset.png", width = 1100, height = 1100)
plot_grid(plot_hl_pre, plot_hl_post, plot_hl_cor_gtex, plot_hl_cor_tcga, labels = LETTERS[1:4], ncol = 2 )
dev.off()


## Low pre - high full gene set - GO:0032274
df_path_lowhi <- data.frame(gtex_pre = weights_gtex_pre["GO:0032274", ], tcga_pre = -weights_tcga_pre["GO:0032274", ], gene = colnames(weights_tcga_pre),
                          gtex = weights_gtex["GO:0032274", ], tcga = weights_tcga["GO:0032274", ]) %>%
              filter(gtex != 0 & tcga != 0) %>%
              mutate(gene = mapIds(org.Hs.eg.db, gene, keytype= "ENSEMBL", column="SYMBOL"))

plot_lh_pre <- df_path_lowhi %>%  ggplot(aes(x = gtex_pre, y = tcga_pre)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_lowhi$gtex_pre, df_path_lowhi$tcga_pre)),
   x = -0.15, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0032274  - step 1") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

genes.lh <- subset(path.map, PathwayID == "GO:0032274")$Symbol
genes.lh_Symbol <- mapIds(org.Hs.eg.db, genes.lh, keytype= "ENSEMBL", column="SYMBOL")

gene_cors.lh_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.lh, ]))))
rownames(gene_cors.lh_gtex ) <- colnames(gene_cors.lh_gtex ) <- genes.lh_Symbol

plot_lh_cor_gtex <-  ggcorrplot(gene_cors.lh_gtex, method = "circle", hc.order = TRUE) +
  ggtitle("GO:0032274 - GTEx") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

gene_cors.lh_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.lh, ]))))
rownames(gene_cors.lh_tcga ) <- colnames(gene_cors.lh_tcga ) <- genes.lh_Symbol

plot_lh_cor_tcga <- ggcorrplot(gene_cors.lh_tcga,  method = "circle", hc.order = TRUE) +
  ggtitle("GO:0032274 - TCGA") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_lh_post <-  df_path_lowhi %>%  ggplot(aes(x = gtex, y = tcga)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data =  data.frame(label = sprintf("r = %.2f", cor(df_path_lowhi$gtex, df_path_lowhi$tcga)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  geom_label_repel(aes(label = gene)) +
  ggtitle("GO:0032274 - full training") +
  xlab("Weights in GTEx") +
  ylab("Weights in TCGA") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

png("figures/data_models_lowhigh_geneset.png", width = 1100, height = 1100)
plot_grid(plot_lh_pre, plot_lh_post, plot_lh_cor_gtex, plot_lh_cor_tcga, labels = LETTERS[1:4], ncol = 2 )
dev.off()
