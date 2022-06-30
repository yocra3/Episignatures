#'#################################################################################
#'#################################################################################
#' Compare
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
library(GOfuncR)
library(GOSim)

## Load data
genes <- read.table("./results/GTEx_coding/input_genes.txt")$V1
path.map <- read.table("results/GTEx_coding/go_kegg_final_gene_map.tsv", header = TRUE)

paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])

weights_gtex <- h5read("results/GTEx_coding/paths_filt2_full_v3.6/model_trained/model_weights.h5","weights_paths")
rownames(weights_gtex) <- paths.vec
colnames(weights_gtex) <- genes

weights_tcga <- h5read("results/TCGA_coding_all/paths_filt2_full_v3.6/model_trained/model_weights.h5","weights_paths")
rownames(weights_tcga) <- paths.vec
colnames(weights_tcga) <- genes

gtex.feat <- read.table("results/TCGA_all/paths_filt2_full_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
tcga.feat <- read.table("results/TCGA_all_TCGA/paths_filt2_full_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
tcga.ctrl.feat <- read.table("results/TCGA_coding_control/paths_filt2_full_v3.6/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(gtex.feat) <- colnames(tcga.feat) <- colnames(tcga.ctrl.feat) <- paths.vec

gtex.feat.pre <- read.table("results/TCGA_all/paths_filt2_pre_v3.8/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
tcga.feat.pre <- read.table("results/TCGA_all_TCGA/paths_filt2_pre_v3.8/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
colnames(gtex.feat.pre) <- colnames(tcga.feat.pre) <- paths.vec

gtex.vst <- loadHDF5SummarizedExperiment("results/GTEx/", prefix = "vst_all_")
tcga.vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

path.N <- table(path.map$PathwayID) %>%
  data.frame()

## Load functions
readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/TCGA_all_TCGA/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
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

path_vals_tcga <-  readPathways("paths_filt2_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)

## Define replicability
tcga.cors <- sapply(paths.vec, pathwayCorr, path_list = path_vals_tcga)
colnames(tcga.cors) <- paths.vec
df.tcga <- makeDFsum(tcga.cors,"full") %>%
  left_join(data.frame(path.N) %>% mutate(path = Var1) %>% dplyr::select(-Var1), by = "path") %>%
  mutate(Database = ifelse(substring(path, 1, 2) == "GO", "GO", "KEGG"))

#

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

path_cors_ctrl <- sapply(seq_len(ncol(gtex.feat)), function(i){

  cor(gtex.feat[, i], tcga.ctrl.feat[, i])
})

path_cors_pre <- sapply(seq_len(ncol(gtex.feat.pre)), function(i){

  cor(gtex.feat.pre[, i], tcga.feat.pre[, i])
})

path_df <- tibble(path = paths.vec, cor_weight = weight_cors, cor_path = path_cors,
                  cor_path_pre = path_cors_pre, cor_path_ctrl = path_cors_ctrl) %>%
  left_join(df.tcga, by = "path")

##
path_df <- mutate(path_df, Category = ifelse(abs(path_cors) > 0.7 & minCor > 0.7 & abs(cor_weight) > 0.7, "High",
        ifelse(abs(path_cors) < 0.3 & minCor > 0.7,
          ifelse(abs(cor_path_pre) < 0.3, "Low pre and full",
              ifelse(abs(cor_path_pre) > 0.7, "Low full", "Other")),
            "Other")))




ggplot(path_df, aes(x = abs(cor_path), y = minCor, color = Category)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))


#
ggplot(path_df, aes(x = abs(cor_path_pre), y = abs(cor_path), color = Category)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))


#
ggplot(path_df, aes(x = Category, y = abs(cor_path), color = Category)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))

ggplot(path_df, aes(x = Category, y = abs(cor_weight), color = Category)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))


#
ggplot(path_df, aes(x = Category, y = abs(cor_path_pre), color = Category)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))


#
ggplot(path_df, aes(x = Category, y = minCor, color = Category)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))

#
ggplot(path_df, aes(x = Category, y = Freq, color = Category)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = c("lightgreen", "lightblue", "blue", "grey"))



cor(abs(path_df$cor_path_pre), abs(path_df$cor_path))


summarize(path_df, n30 = mean(abs(cor_path) > 0.3),
                  n50 = mean(abs(cor_path) > 0.5),
                  n70 = mean(abs(cor_path) > 0.7))

#
summarize(path_df, n30 = mean(abs(cor_path_ctrl) > 0.3),
                  n50 = mean(abs(cor_path_ctrl) > 0.5),
                  n70 = mean(abs(cor_path_ctrl) > 0.7))



summarize(path_df, n30 = mean(abs(cor_path_pre) > 0.3),
                  n50 = mean(abs(cor_path_pre) > 0.5),
                  n70 = mean(abs(cor_path_pre) > 0.7))

summarize(path_df, n30 = mean(abs(minCor) > 0.3),
                  n50 = mean(abs(minCor) > 0.5),
                  n70 = mean(abs(minCor) > 0.7))



# ggplot(path_df, aes(x = Freq, y = abs(cor_path))) +
#   geom_point() +
#   theme_bw()
# cor(abs(path_df$cor_path), path_df$Freq)
#
ggplot(path_df, aes(x = abs(path_cors), y = minCor, color = abs(cor_weight))) +
  geom_point() +
  theme_bw()


# ggplot(path_df, aes(x = Category, y = Freq)) +
#   geom_boxplot() +
#   theme_bw()

## Correlation in high correlation pathways
genes.top <- subset(path.map, PathwayID == "GO:0003356")$Symbol
gene_cors.top_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.top, ]))))
corrplot(gene_cors.top_gtex, method = "number", order = "hclust")

gene_cors.top_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.top, ]))))
corrplot(gene_cors.top_tcga, method = "number", order = "hclust")


genes.low <- subset(path.map, PathwayID == "GO:0030575")$Symbol
gene_cors.low_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.low, ]))))
corrplot(gene_cors.low_gtex, method = "number", order = "hclust")

gene_cors.low_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.low, ]))))
corrplot(gene_cors.low_tcga, method = "number", order = "hclust")

genes.low2 <- subset(path.map, PathwayID == "GO:0035751")$Symbol
gene_cors.low2_gtex <- cor(t(data.matrix(assay(gtex.vst[genes.low2, ]))))
corrplot(gene_cors.low2_gtex, method = "number", order = "hclust")

gene_cors.low2_tcga <- cor(t(data.matrix(assay(tcga.vst[genes.low2, ]))))
corrplot(gene_cors.low2_tcga, method = "number", order = "hclust")



topPaths <- subset(path_df, Category == "High")$path
badPaths1 <- subset(path_df, Category == "Low pre and full")$path
badPaths2 <- subset(path_df, Category == "Low full")$path

#
# ## Compute GO terms similarities
# sims_mat <- getTermSim(paths.vec, method = "Lin", verbose = TRUE)
# sims_mat_filt <- sims_mat[rowMeans(sims_mat == 0) < 1, colMeans(sims_mat == 0) < 1]
#
# sim_top <- sims_mat_filt[topPaths[topPaths %in%  rownames(sims_mat_filt)], topPaths[topPaths %in%  rownames(sims_mat_filt)]]
# sim_low <- sims_mat_filt[badPaths[badPaths %in%  rownames(sims_mat_filt)], badPaths[badPaths %in%  rownames(sims_mat_filt)]]
# sim_top_low <- sims_mat_filt[topPaths[topPaths %in%  rownames(sims_mat_filt)], badPaths[badPaths %in%  rownames(sims_mat_filt)]]
# sim_rest <- sims_mat_filt[!rownames(sims_mat_filt) %in% c(topPaths,badPaths), !rownames(sims_mat_filt) %in% c(topPaths,badPaths)]
#
#
#
# sim_df <- data.frame(scores = c(sim_top[upper.tri(sim_top)], sim_low[upper.tri(sim_low)], as.vector(sim_top_low), sim_rest[upper.tri(sim_rest)]),
#                       Type = rep(c("Top", "Low", "Top-Low", "All"), lengths(list(sim_top[upper.tri(sim_top)], sim_low[upper.tri(sim_low)], as.vector(sim_top_low), sim_rest[upper.tri(sim_rest)])))) %>%
#                       mutate(scores = ifelse(scores > 1, 1, scores))
#
# ggplot(sim_df, aes(x = Type, y = scores)) + geom_boxplot()
#
#
# sim_both <- 1 - sims_mat_filt[ rownames(sims_mat_filt) %in% c(topPaths,badPaths) , rownames(sims_mat_filt) %in% c(topPaths,badPaths)]
# sim_both[sim_both < 0 ] <- 0
# sim_clust <- hclust(as.dist(sim_both))
#
# library(dendextend)
# dend <- as.dendrogram(sim_clust)
#
# col <- ifelse(colnames(sim_both) %in% topPaths, "red", "blue")
#
# labels_colors(dend) <- col[order.dendrogram(dend)]
# plot(dend)
#
# groups <- cutree(sim_clust, h = 0.99)
# json_l <- lapply(unique(groups), function(i){
#   l <- lapply(colnames(sim_both)[groups == i], function(x){
#     if (x %in% topPaths){
#       list(fill = "#fcba03")
#     }  else if (x %in% badPaths){
#           list(fill = "#03bafc")
#         }
#   })
#   names(l) <- colnames(sim_both)[groups == i]
#   lj <- toJSON(l)
#   write(lj, paste0("results/GTEx_coding/cluster", i, ".json") )
# })
#
# outPaths <- c(topPaths, badPaths)
# outPaths <- outPaths[!outPaths %in% colnames(sim_both)]
# l2 <- lapply(outPaths, function(x){
#   if (x %in% topPaths){
#     list(fill = "#fcba03")
#   }  else if (x %in% badPaths){
#         list(fill = "#03bafc")
#       }
# })
# names(l2) <- outPaths
# lj <- toJSON(l2)
# write(lj, "results/GTEx_coding/outside_cluster.json")


all_gos <- get_child_nodes("GO:0008150")
top_gos <- c(subset(all_gos, distance < 3)$child_go_id, "GO:0006810", "GO:0048731", "GO:0048513",
  "GO:0006725", "GO:0016043", "GO:0071310", "GO:0010033", "GO:0010646", "GO:0071705", "GO:0043603", "GO:0007399",
  "GO:0043207", "GO:1901135", "GO:0009888","GO:0035295", "GO:0044249", "GO:0044248", "GO:0051173", "GO:0030154","GO:0051128", "GO:0008104",
   "GO:0031323", "GO:0032535", "GO:0060284","GO:0048638", "GO:0030030", "GO:0006811", "GO:0010817"     )


selPaths <-  c(topPaths, badPaths1, badPaths2)
selPaths <- selPaths[grep("GO", selPaths)]
go_clust <- list(selPaths[1])
selPaths <- selPaths[-1]
for (path in selPaths){
  used <- FALSE
  for (i in seq_len(length(go_clust))){
    anc <- getMinimumSubsumer(path, go_clust[[i]][[1]])
    if (! anc %in% top_gos) {
      go_clust[[i]] <- c(go_clust[[i]], path)
      used <- TRUE
    }
  }
  if(!used){
    go_clust <- c(go_clust, path)
  }
}
go_clust_filt <- go_clust[lengths(go_clust) > 2]

json_l <- lapply(seq_len(length(go_clust_filt)), function(i){
  l <- lapply(go_clust_filt[[i]], function(x){
    if (x %in% topPaths){
      list(fill = "#fcba03")
    }  else if (x %in% badPaths1){
          list(fill = "#03bafc")
        }

        else if (x %in% badPaths2){
              list(fill = "#03f4fc")
            }
    else {
      list(fill = "#bbbbbb")
    }
  })
  names(l) <- go_clust_filt[[i]]
  lj <- toJSON(l)
  write(lj, paste0("results/GTEx_coding/cluster", i, ".json") )
})



rem_paths <- c("GO:0060603", "GO:0090189", "GO:2001212", "GO:0022012", "GO:0021846", "GO:0021534",
"GO:0021952", "GO:0034390", "GO:0097050", "GO:0097048", "GO:0045019", "GO:1903299", "GO:0006085",
"GO:1901881", "GO:0010592", "GO:0045624", "GO:0045625", "GO:0048843", "GO:0014829",
"GO:0071472", "GO:0036295", "GO:0110096", "GO:0061481", "GO:0043619", "GO:0002638", "GO:0010955",
"GO:0035278", "GO:1903800", "GO:1900152", "GO:0031445", "GO:0140719", "GO:0071514", "GO:1902915",
"GO:0033234", "GO:1901984", "GO:0031061", "GO:0050686", "GO:0042532", "GO:1902902", "GO:0043508", "GO:0006123",
"GO:0031998", "GO:0070125", "GO:0031145", "GO:0006515", "GO:0042159", "GO:1904292", "GO:0071712", "GO:0030214",
"GO:0007039", "GO:0035815", "GO:2000193", "GO:2000833", "GO:2001138", "GO:0006415", "GO:0032986", "GO:1904896", "GO:0010715",
"GO:0045050", "GO:0045048", "GO:0042407", "GO:0051561", "GO:0099566", "GO:0003091", "GO:0001780", "GO:0030277",
"GO:0031649", "GO:0061484", "GO:0033561", "GO:1901386", "GO:0003096", "GO:0043201", "GO:1902065",
"GO:0071389", "GO:0038166", "GO:0007191", "GO:0060159", "GO:0034497", "GO:0099638", "GO:0007144", "GO:0051307",
"GO:0046599", "GO:0035721", "GO0048753", "GO:0060155", "GO:0033127", "GO:0035404", "GO:0001516", "GO:0045540", "GO:0019372",
"GO:0009263", "GO:0030810", "GO:0090208", "GO:0010866", "GO:0046339", "GO:0090153", "GO:0031282", "GO:2001169",
"GO:2001267", "GO:1905048", "GO:0051798", "GO:0048820",  "GO:0032525", "GO:0060670", "GO:1900246", "GO:0090189",  "GO:2001212", "GO:0061626", "GO:0060977",
"GO:0021542","GO:0021548", "GO:0097050", "GO:0097048", "GO:1902337", "GO:0034390", "GO:0044346", "GO:1901028", "GO:1901030", "GO:0043653",
"GO:0090201", "GO:1903299", "GO:0045019","GO:0006085", "GO:0045624", "GO:0045625", "GO:0014829", "GO:1902855", "GO:1902018", "GO:0045724", "GO:0010592",
"GO:0070584", "GO:0070170", "GO:0045605", "GO:0003334", "GO:0045606", "GO:0060221", "GO:0042670", "GO:0042481", "GO:0002092",
"GO:0033089", "GO:0045579", "GO:0035485", "GO:1903975", "GO:0072673", "GO:0030852", "GO:0045649", "GO:0045655",
"GO:0061052", "GO:0045603", "GO:2000696", "GO:0070307", "GO:0001886", "GO:0090557", "GO:0060117", "GO:0032332", "GO:0061037", "GO:0002068",
"GO:0010831", "GO:0010832", "GO:0002070", "GO:0097084", "GO:0010635", "GO:0045793", "GO:1905208", "GO:0002052",
"GO:1901881", "GO:0090141", "GO:0045672", "GO:0048843", "GO:0070593", "GO:0071679", "GO:0036035", "GO:2001198", "GO:0022011",
"GO:2000737","GO:2000738", "GO:0045161", "GO:0016322", "GO:1902285", "GO:0051152", "GO:0045654","GO:2001223","GO:2001224", "GO:2001014",
"GO:0045662", "GO:0061003", "GO:0070572", "GO:0035855", "GO:0061000", "GO:0061418", "GO:0043619", "GO:0036295",
"GO:0060397", "GO:0035747", "GO:2000508", "GO:0060394", "GO:0030949", "GO:0030948", "GO:0038065", "GO:0023035", "GO:0038007",
"GO:0071472", "GO:1900037", "GO:1903209", "GO:0002689", "GO:0035767", "GO:0061314", "GO:0003306", "GO:0060391",
"GO:0043567", "GO:0038092", "GO:0032925", "GO:0032926", "GO:0110096", "GO:0035728", "GO:0036119")
selPaths <- paths.vec[!paths.vec %in% c(topPaths, badPaths,rem_paths)]
selPaths <- selPaths[grep("GO", selPaths)]

for (path in selPaths){
  for (i in seq_len(length(go_clust_filt))){
    anc <- getMinimumSubsumer(path, go_clust_filt[[i]][[1]])
    if (! anc %in%   top_gos      ){
      go_clust_filt[[i]] <- c(go_clust_filt[[i]], path)
      break
    }
  }
}








getMinimumSubsumer(path, go_clust_filt[[i]][[1]])
