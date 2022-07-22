#'#################################################################################
#'#################################################################################
#' Select GO and KEGG features for model
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(tidyverse)
library(parallel)
library(matrixStats)
library(ggrepel)
library(ggcorrplot)
library(cowplot)
# library(GOSemSim)

## Load functions
readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    tab <- read.table(path, header = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}

pathwayCorr <- function(path_list, col){
  path_mat <- sapply(path_list, function(x) x[, col])
  cors <- cor(path_mat)
  cors[upper.tri(cors)]
}


makeDFsum <- function(cors, mod_name){
  df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
    mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
            model = mod_name) %>%
            as_tibble()

}


## Load pathways annotations
kegg.map <- read.table("results/preprocess/go_kegg_gene_map.tsv", header = TRUE)
paths <- read.table("results/GTEx_coding/comb_paths_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
input_genes <- read.table("results/GTEx_coding/input_genes.txt", header = FALSE)



kegg.map.com <- subset(kegg.map, PathwayID %in% paths.vec & Symbol %in% input_genes$V1)
kegg.N <- table(kegg.map.com$PathwayID) %>%
  data.frame()

## Define similarity between pathways
path_genes <- mclapply(paths.vec, function(x) subset(kegg.map.com, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
# mean_l <- mclapply(paths.vec, function(x) sapply(paths.vec, function(y) mean(path_genes[[x]] %in% path_genes[[y]])), mc.cores = 20)
# mean_mat <- matrix(unlist(mean_l), nrow = length(paths.vec))
# rownames(mean_mat) <- colnames(mean_mat) <- paths.vec

gene_mat <- matrix(0, length(paths.vec), length(unique( kegg.map.com$Symbol)), dimnames = list(paths.vec, unique(kegg.map.com$Symbol)))
 for (i in paths.vec) gene_mat[i, path_genes[[i]]] <- 1
gene_d <- dist(gene_mat, "binary")
save(gene_d, file = "results/GTEx_coding/go_kegg_pathways_distance.Rdata")
gene_dmat <- as.matrix(gene_d)

# hsGO <- godata('org.Hs.eg.db', ont="BP")
# sim_go <- mgoSim(good_paths, good_paths, semData=hsGO, measure="Wang", combine=NULL)

## Load values of pathways without full training
path_vals_full <- readPathways("comb_paths_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
path_vals_pre <- readPathways("comb_paths_v3.8", sufix = c("", letters[1:5]), path_name = paths.vec)

## Define replicability
full.cors <- sapply(paths.vec, pathwayCorr, path_list = path_vals_full)
colnames(full.cors) <- paths.vec
df.full <- makeDFsum(full.cors,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))


pre.cors <- sapply(paths.vec, pathwayCorr, path_list = path_vals_pre)
colnames(pre.cors) <- paths.vec
df.pre <- makeDFsum(pre.cors,"pretrained") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

png("figures/N_vs_minCor_pathways.png")
ggplot(df.pre, aes(x = Freq, y = minCor, color = Database)) +
  geom_point() +
  scale_x_log10(name = "Genes per pathway") +
  scale_y_continuous(name = "Replicability") +
  theme_bw()
dev.off()

png("figures/N_vs_minCor_pathways_categorical.png")
df.pre %>%
  mutate(Freq_groups = cut(Freq, c(0, 20, 50, 100, 200, 300, 1000, 100000), c("<20", "20-50", "50-100", "100-200", "200-300", "300-1000", ">1000"))) %>%
  ggplot(aes(x = Freq_groups, y = minCor)) +
    geom_boxplot() +
    scale_x_discrete(name = "Genes per pathway") +
    scale_y_continuous(name = "Replicability") +
    theme_bw()
dev.off()


png("figures/N_vs_minCor_pathways_fulltrained.png")
ggplot(df.full, aes(x = Freq, y = minCor, color = Database)) +
  geom_point() +
  scale_x_log10(name = "Genes per pathway") +
  scale_y_continuous(name = "Replicability") +
  theme_bw()
dev.off()

png("figures/N_vs_minCor_pathways_categorical_fulltrained.png")
df.full %>%
  mutate(Freq_groups = cut(Freq, c(0, 20, 40, 60, 80, 100, 300, 100000), c("<20", "20-40", "40-60", "60-80", "80-100", "100-300",  ">300"))) %>%
  ggplot(aes(x = Freq_groups, y = minCor)) +
    geom_boxplot() +
    scale_x_discrete(name = "Genes per pathway") +
    scale_y_continuous(name = "Replicability") +
    theme_bw()
dev.off()

df.full2 %>%
  mutate(Freq_groups = cut(Freq, c(0, 20, 30, 40, 60), c("<20", "20-30", "30-40", "40+"))) %>%
  ggplot(aes(x = Freq_groups, y = minCor)) +
    geom_boxplot() +
    scale_x_discrete(name = "Genes per pathway") +
    scale_y_continuous(name = "Replicability") +
    theme_bw()

df.comb <- rbind(df.pre, df.full)  %>%
  dplyr::select(-medCor, -class) %>%
  spread(model, minCor)

png("figures/minCor_pathways_full_vs_pretrained.png")
df.comb %>%
  mutate(group = ifelse(Freq > 100, "> 100 genes", "< 100 genes")) %>%
  ggplot(aes(x = pretrained, y = full, col = group)) +
   geom_point() +
   theme_bw() +
   scale_x_continuous(name = "Replicability pretrained model") +
   scale_y_continuous(name = "Replicability full model")
dev.off()

## Select good pathways by correlation
good_paths <- as.character(subset(kegg.N, Freq < 30 & Freq >= 10)$Var)
## Compute correlation between pathways
cor_l <- mclapply(path_vals_pre, function(x) cor(x[, good_paths]), mc.cores = 6)
cor_mat <- sapply(cor_l, function(x) x[upper.tri(x)])
df.cor <- tibble(minCor = rowMins(abs(cor_mat)), meanCor = rowMeans(abs(cor_mat)),
  dist = gene_dmat_sel[upper.tri(gene_dmat_sel)])

cor_mat0 <- cor_l[[1]]
rownames(cor_mat0) <- colnames(cor_mat0) <- good_paths
diag(cor_mat0) <- 0
cor_df <- tibble(path = rownames(cor_mat0), max_cor = rowMaxs(abs(cor_mat0)))


cor_full0 <- cor(path_vals_full[[1]][, good_paths])

gene_dmat_sel <- gene_dmat[good_paths, good_paths]
gene_dmat_sel2 <- gene_dmat_sel
diag(gene_dmat_sel2) <- 1
dist_df <- tibble(path = rownames(gene_dmat_sel2), min_dist = rowMins(gene_dmat_sel2))

df.comb.sel <- filter(df.comb, path %in% good_paths) %>%
  mutate(diff_rep = pretrained - full) %>%
  left_join(dist_df, by = "path") %>%
  left_join(sim_df, by = "path") %>%
  left_join(cor_df, by = "path")


png("figures/cor_vs_dist_pathways.png")
df.cor %>%
  mutate(group = ifelse(dist < 1, "Shared Genes", "No common genes")) %>%
  ggplot(aes(x = group, y = meanCor)) +
    geom_boxplot() +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Pathways correlation (absolute value)")
dev.off()


png("figures/cor_vs_dist_pathways_sharedGenes.png")
df.cor %>%
  filter(dist < 1) %>%
  mutate(cat = cut(dist, seq(0, 1, 0.1), labels = paste0(seq(0, 90, 10), "-", seq(10, 100, 10), "%"), include.lowest = TRUE )) %>%
  ggplot(aes(x = cat, y = meanCor)) +
    geom_boxplot() +
    theme_bw() +
    scale_x_discrete(name = "Proportion of unique genes") +
    scale_y_continuous(name = "Pathways correlation (absolute value)")
dev.off()

summary(lm(meanCor ~ dist, df.cor))
summary(lm(meanCor ~ log10(dist + 1e-3), df.cor))

summary(lm(meanCor ~ dist, df.cor, subset = dist < 1))

summary(lm(meanCor ~ log10(dist + 1e-3), df.cor, subset = dist < 1))

gene_dmini <- as.dist(gene_dmat_sel)
cor_mean <- matrix(1, length(good_paths), length(good_paths))
cor_mean[upper.tri(cor_mean)] <- df.cor$meanCor
cor_mean[lower.tri(cor_mean)] <- t(cor_mean)[lower.tri(cor_mean)]


colnames(cor_mean) <- rownames(cor_mean) <- colnames(gene_dmat_sel)

cuts <- seq(0.1, 0.9, 0.1)
minCors_list <- lapply(cuts, function(y){
  cls <- cutree(hclust(gene_dmini), h = y)
  minCors <- sapply(unique(cls), function(x) {
    minicor <- cor_mean[cls == x, cls == x]
    if (!is(minicor, "matrix")){
      return(NA)
    } else {
      min(abs(minicor))
    }
  })
  minCors
})
minCors_df <- data.frame(minCor = unlist(minCors_list),
                          cut = rep(cuts, lengths(minCors_list)))

png("figures/clustering_distance_vs_cor.png")
minCors_df %>% ggplot(aes(x = factor(cut), y = minCor)) + geom_boxplot() + theme_bw()
dev.off()

minCors_df %>%
  group_by(cut) %>%
  summarize(n = n(), nclust = sum(!is.na(minCor)), mean90 = mean(minCor > 0.9, na.rm = TRUE),
            mean80 = mean(minCor > 0.8, na.rm = TRUE),
            mean70 = mean(minCor > 0.7, na.rm = TRUE))

## Cluster pathways until all have a binary distance < 0.5
sel_paths <- good_paths
gene_dloop <- gene_dmini
gene_dmat_loop <- gene_dmat_sel
path_cls <- cutree(hclust(gene_dmini), h = 0.5)
while(length(sel_paths) != length(unique(path_cls))){
  sel_paths <- sapply(unique(path_cls), function(cl){
    paths <- rownames( gene_dmat_loop)[path_cls == cl]
    df.sub <- subset(kegg.N, Var1 %in% paths)
    as.character(df.sub$Var1[which.max(df.sub$Freq)])
  })
  gene_dmat_loop <- gene_dmat_loop[sel_paths, sel_paths]
  gene_dloop <- as.dist(gene_dmat_loop)
  path_cls <- cutree(hclust(gene_dloop), h = 0.5)
}
gene_dmat_filt <- gene_dmat[sel_paths, sel_paths]
diag(gene_dmat_filt) <- 1
summary(rowMins(gene_dmat_filt))

kegg.map.filt <- subset(kegg.map, PathwayID %in% sel_paths)
write.table(kegg.map.filt, file = "results/GTEx_coding/go_kegg_filt_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


##
paths2 <- read.table("results/GTEx_coding/paths_filt1_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec2 <- as.character(paths2[, 1])

path_vals_full2 <- readPathways("paths_filt1_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec2)

full.cors2 <- sapply(paths.vec2, pathwayCorr, path_list = path_vals_full2)
colnames(full.cors2) <- paths.vec2
df.full2 <- makeDFsum(full.cors2,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))



path_vals_full2b <- readPathways("paths_filt1_full_v3.11", sufix = c("", letters[1:5]), path_name = paths.vec2)

full.cors2b <- sapply(paths.vec2, pathwayCorr, path_list = path_vals_full2b)
colnames(full.cors2b) <- paths.vec2
df.full2b <- makeDFsum(full.cors2b,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))


sel_paths2 <- subset(df.full2b, minCor > 0.7)$path
kegg.map.filt2 <- subset(kegg.map, PathwayID %in% sel_paths2 & Symbol %in% input_genes$V1)
write.table(kegg.map.filt2, file = "results/GTEx_coding/go_kegg_filt2_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


##
paths3 <- read.table("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec3 <- as.character(paths3[, 1])

path_vals_full3 <- readPathways("paths_filt2_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec3)


full.cors3 <- sapply(paths.vec3, pathwayCorr, path_list = path_vals_full3)
colnames(full.cors3) <- paths.vec3
df.full3 <- makeDFsum(full.cors3,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

#
paths3b <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec3b <- as.character(paths3b[, 1])

path_vals_full3b <- readPathways("paths_filt2_full_v3.11", sufix = c("", letters[1:5]), path_name = paths.vec3b)


full.cors3b <- sapply(paths.vec3b, pathwayCorr, path_list = path_vals_full3b)
colnames(full.cors3b) <- paths.vec3b
df.full3b <- makeDFsum(full.cors3b,"full") %>%
  left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

full.cors3b_main <- sapply(paths.vec3b, function(x) median(abs(cor(sapply(path_vals_full3b, function(y) y[, x]))[1,-1])))

df.3b_worse <- df.full3b %>%
  mutate(main_cor = full.cors3b_main) %>%
  filter(minCor < 0.7)

png("figures/select_features_worse.png")
ggplot(df.3b_worse, aes( x = minCor, y = main_cor) ) + geom_point() +
  theme_bw() +
  xlab("Replicability") +
  ylab("Median correlation main model") +
  geom_label_repel(data = subset(df.3b_worse, main_cor < 0.7), aes(label = path))
dev.off()

#
# sel_paths3 <- subset(df.full3, minCor > 0.7)$path
# kegg.map.filt3 <- subset(kegg.map, PathwayID %in% sel_paths3 & Symbol %in% input_genes$V1)
# write.table(kegg.map.filt3, file = "results/GTEx_coding/go_kegg_filt3_gene_map.tsv",
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

mat1 <- sapply(path_vals_full3b, function(y) y[, "path:hsa00591"])
cor1 <- cor(mat1)
rownames(cor1) <- colnames(cor1) <- c("Main", paste("Initialization", 1:5))
plot_cor1 <- ggcorrplot(cor1,  method = "circle") +
  ggtitle("path:hsa00591") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_path1 <- data.frame(Main = mat1[, 1], mod = mat1[, 4]) %>%
  ggplot(aes( x = Main, y = mod)) +
  geom_point() +
  theme_bw() +
  ggtitle("path:hsa00591") +
  xlab("Initialization 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))


png("figures/select_features_hsa00591.png", width = 1100, height = 500)
plot_grid(plot_cor1, plot_path1, labels = LETTERS[1:2], ncol = 2 )
dev.off()

mat2 <- sapply(path_vals_full3b, function(y) y[, "GO:0045663"])
cor2 <- cor(mat2)
rownames(cor2) <- colnames(cor2) <- c("Main", paste("Initialization", 1:5))
plot_cor2 <- ggcorrplot(cor2,  method = "circle") +
  ggtitle("GO:0045663") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_path2 <- data.frame(Main = mat2[, 1], mod = mat2[, 4]) %>%
  ggplot(aes( x = Main, y = mod)) +
  geom_point() +
  theme_bw() +
  ggtitle("GO:0045663") +
  xlab("Initialization 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

png("figures/select_features_GO0045663.png", width = 1100, height = 500)
plot_grid(plot_cor2, plot_path2, labels = LETTERS[1:2], ncol = 2 )
dev.off()

mat3 <- sapply(path_vals_full3b, function(y) y[, "GO:0042538"])
cor3 <- cor(mat3)
rownames(cor3) <- colnames(cor3) <- c("Main", paste("Initialization", 1:5))
plot_cor3 <- ggcorrplot(cor3,  method = "circle") +
  ggtitle("GO:0042538") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

plot_path3 <- data.frame(Main = mat3[, 1], mod = mat3[, 2]) %>%
  ggplot(aes( x = Main, y = mod)) +
  geom_point() +
  theme_bw() +
  ggtitle("GO:0042538") +
  xlab("Initialization 1") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

png("figures/select_features_GO0042538.png", width = 1100, height = 500)
plot_grid(plot_cor3, plot_path3, labels = LETTERS[1:2], ncol = 2 )
dev.off()

#
# #
# paths4 <- read.table("results/GTEx_coding/paths_filt3_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
# paths.vec4 <- as.character(paths4[, 1])
#
# path_vals_full4 <- readPathways("paths_filt3_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec4)
#
#
# full.cors4 <- sapply(paths.vec4, pathwayCorr, path_list = path_vals_full4)
# colnames(full.cors4) <- paths.vec4
# df.full4 <- makeDFsum(full.cors4, "full") %>%
#   left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
#   mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))
