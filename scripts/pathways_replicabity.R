#'#################################################################################
#'#################################################################################
#' Explore replicability in pathways activations from different models
#'#################################################################################
#'#################################################################################


## Load libraries ####
# library(topGO)
library(SummarizedExperiment)
library(tidyverse)
library(rjson)
library(rhdf5)
library(HDF5Array)
library(rtracklayer)



readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    # path <- paste0("results/TCGA_gexp_combat_std/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    path <- paste0("results/TCGA_gexp_combat_coding_std/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    tab <- read.table(path, header = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}


readPathways2 <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    # path <- paste0("results/TCGA_gexp_combat_std/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    path <- paste0("results/TCGA_gexp_combat_coding_std/", model, i, "/model_features/prune_low_magnitude_dense_1.tsv")
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

cancer.keg <- c("hsa04010", "hsa04310", "hsa04350", "hsa04370", "hsa04630", "hsa04024", "hsa04151", "hsa04150", "hsa04110", "hsa04210", "hsa04115", "hsa04510", "hsa04520", "hsa03320")

kegg.map <- read.table("results/preprocess/kegg_filt_gene_map.tsv", header = TRUE)
kegg.N <- table(kegg.map$PathwayID)

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

# paths <- read.table("results/TCGA_gexp_norm/kegg_filt_v3.2/model_trained/pathways_names.txt", header = TRUE)
paths <- read.table("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2/model_trained/pathways_names.txt", header = TRUE)

paths.vec <- as.character(paths[, 1])
kegg.df.com <- subset(kegg.df, pathID %in% paths.vec)
kegg.genes.N <- kegg.map %>%
  group_by(Symbol) %>%
  summarize(N = n())

## Load correlations
readCors <- function(base, models, paths.name, model.name){
  base <- readPathways(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}
base <- readCors("kegg_filt2_v3.5", c("", letters[1:5]), paths.vec, "Pathway") %>%
  mutate(training = "primed")
base2 <- readCors("kegg_filt2_v3.6", c("", letters[1:5]), paths.vec, "Pathway") %>%
  mutate(training = "primed + pretrained")
base3 <- readCors("kegg_filt2_v3.7", c("", letters[1:5]), paths.vec, "Pathway") %>%
  mutate(training = "primed + dropout")

post <- readCors("kegg_filt2_v4.2", c("", letters[1:5]), paths.vec, "Pathway + Dense")  %>%
  mutate(training = "primed")
post2 <- readCors("kegg_filt2_v4.3", c("", letters[1:5]), paths.vec, "Pathway + Dense")  %>%
  mutate(training = "primed + pretrained")
post3 <- readCors("kegg_filt2_v4.4", c("", letters[1:5]), paths.vec, "Pathway + Dense")  %>%
  mutate(training = "primed + dropout")


readCors2 <- function(base, models, paths.name, model.name){
  base <- readPathways2(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}
pre <- readCors2("kegg_filt2_v6.1", c("", letters[1:5]), paths.vec, "Dense + Pathway") %>%
  mutate(training = "primed")
pre2 <- readCors2("kegg_filt2_v6.2",c("", letters[1:5]), paths.vec, "Dense + Pathway") %>%
  mutate(training = "primed + pretrained")
pre3 <- readCors2("kegg_filt2_v6.3", c("", letters[1:5]), paths.vec, "Dense + Pathway") %>%
  mutate(training = "primed + dropout")

pre.post <- readCors2("kegg_filt2_v5.2", c("", letters[1:5]), paths.vec, "Dense + Pathway + Dense") %>%
  mutate(training = "primed")
pre2.post <- readCors2("kegg_filt2_v5.3", c("", letters[1:5]), paths.vec, "Dense + Pathway + Dense") %>%
  mutate(training = "primed + pretrained")
pre3.post <- readCors2("kegg_filt2_v5.4", c("", letters[1:5]), paths.vec, "Dense + Pathway + Dense") %>%
  mutate(training = "primed + dropout")

df.all <- Reduce(rbind, list(base, base2, base3, post, post2, post3, pre, pre2, pre3, pre.post, pre2.post , pre3.post )) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(model = factor(model , levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")),
          training = factor(training , levels = c("primed", "primed + pretrained", "primed + dropout")))



# path_vals <- readPathways("kegg_filt2_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
# path_vals_pre <- readPathways2("kegg_filt2_v6.2", sufix = c("", letters[1:5]), path_name = paths.vec)
#
# sel.cors <- sapply(paths.vec, pathwayCorr, path_list = c(path_vals, path_vals_pre))

## Plots
png("figures/minCor_models.png", width = 1500)
df.all %>%
  ggplot(aes(x = model, y = minCor, color = model)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ training) +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
dev.off()


png("figures/minCor_pathwayCat.png",width = 1000)
df.all %>%
  filter(training ==  "primed + pretrained")   %>%
  filter(!is.na(top_cat)) %>%
  ggplot(aes(x = top_cat, y = minCor, color = model)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
dev.off()

png("figures/minCor_pathwayCat2.png", width = 3000, height = 1000)
df.all %>%
  filter(training ==  "primed + pretrained")   %>%
  filter(!is.na(category)) %>%
  filter(top_cat != "Environmental Information Processing") %>%
  ggplot(aes(x = category, y = minCor, color = model)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~top_cat, scales = "free_x") +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
dev.off()

png("figures/minCor_vs_pathwayNgenes.png", width = 800, height = 800)
df.all %>%
  filter(training ==  "primed + pretrained" & !is.na(top_cat))   %>%
  ggplot(aes(y = minCor, x = Freq, color = top_cat)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~model) +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
dev.off()


png("figures/medCor_models.png", width = 500)
df.all %>%
  ggplot(aes(x = model, y = medCor, color = model)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ training)
dev.off()

#'#################################################################################
## Dropouts
#'#################################################################################
## Get dfs for all dropouts
# cors.drops <- lapply(10:13, function(x){
#   mod <- paste0("kegg_filt_v2.", x)
#   drop20 <- readPathways(mod, letters[1:5], paths.vec)
#
#   cors20 <- sapply(paths.vec, pathwayCorr, path_list = drop20)
#   colnames(cors20) <- paths.vec
#   rownames(cors20 ) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
#
#   cors20
# })
#
# ## Get dfs for all dropouts
# df.drops <- Map(makeDFsum, cors.drops, paste("Dropout", c("20%", "50%", "75%", "90%")))
#
# df.drops$cot <- df.base
# df.all <- Reduce(rbind, df.drops ) %>%
#   left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path")
#
# ## Plots
# png("figures/minCor_models.png", width = 500)
# df.all %>%
#   ggplot(aes(x = model, y = minCor, color = model)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
#
#
# png("figures/minCor_pathwayCat.png", width = 1000)
# df.all %>%
#   ggplot(aes(x = top_cat, y = minCor, color = model)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
#
# png("figures/minCor_pathwayCat2.png", width = 3000, height = 1000)
# df.all %>%
#   ggplot(aes(x = category, y = minCor, color = model)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_wrap(~top_cat, scales = "free_x", nrow = 3)
# dev.off()
#
# png("figures/medCor_models.png", width = 500)
# df.all %>%
#   ggplot(aes(x = model, y = medCor, color = model)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
#
#
# png("figures/medCor_pathwayCat.png", width = 1000)
# df.all %>%
#   ggplot(aes(x = top_cat, y = medCor, color = model)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
#
# png("figures/medCor_pathwayCat2.png", width = 3000, height = 1000)
# df.all %>%
#   ggplot(aes(x = category, y = medCor, color = model)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_wrap(~top_cat, scales = "free_x", nrow = 3)
# dev.off()
#
#
# summary(lm(minCor ~ model, df.all))
# summary(lm(medCor ~ model, df.all))
#
# df.sum <- df.all %>% group_by(path) %>% summarize(nlow = sum(class == "low"), nhigh = sum(class == "High"))
# table(low = df.sum$nlow, high = df.sum$nhigh)
# badPathways <- subset(df.sum, nlow > 0)$path
# subset(kegg.df.com, pathID %in% badPathways)

## Compute mse
ori <- h5read("results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5","methy")
base <- h5read("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_features/autoencoder_output.h5", "auto")
post <- h5read("results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_features/autoencoder_output.h5", "auto")
pre <- h5read("results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_features/autoencoder_output.h5", "auto")
pre_post <- h5read("results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_features/autoencoder_output.h5", "auto")
mat_list <- list(base, post, pre, pre_post)
names(mat_list) <- c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")
vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_combat_coding/", prefix = "vsd_norm")

labels <- as.character(read.table("./results/TCGA_gexp_combat_coding/individuals_labels.txt")$V1)
genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
rownames(ori) <- rownames(base) <- rownames(post) <- rownames(pre) <- rownames(pre_post) <- as.character(genes$V1)

phenos <- TCGAbiolinks:::get_IDs(vst)


kegg.map.com <- kegg.map %>%
  left_join(mutate(kegg.df.com, PathwayID = pathID) %>%
            select(PathwayID, top_cat), by = "PathwayID")

kegg.genes.cats <- kegg.map.com %>%
  group_by(Symbol) %>%
  summarize(Metabolism = sum(top_cat == "Metabolism"),
            Cellular = sum(top_cat == "Cellular Processes"),
            Environment = sum(top_cat == "Environmental Information Processing"),
            Genetic = sum(top_cat == "Genetic Information Processing"),
            Disease = sum(top_cat == "Human Diseases"),
            Organism = sum(top_cat == "Organismal Systems")
          )


# annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
# annot <- filter(annot, type == "gene") %>%
#   as_tibble() %>%
#      mutate(Symbol = gsub("\\.[0-9]*_[0-9]*", "", gene_id , perl = TRUE))
## Correlation of genes not in model vs genes in model
gens.cor <- cor(t(ori[rownames(ori) %in%  kegg.map$Symbol,]), t(ori[!rownames(ori) %in%  kegg.map$Symbol,]))
summary(colMaxs(abs(gens.cor)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04249 0.53749 0.63315 0.62138 0.72692 0.96761

mean.mses <- sapply(mat_list, function(x) mean((x - ori)**2))
median.mses <- sapply(mat_list, function(x) median((x - ori)**2))
genes.mse <-  sapply(mat_list, function(x) rowMeans((x - ori)**2))
genes.cors <- sapply(mat_list, function(x) sapply(seq_len(nrow(ori)), function(i) cor(ori[i, ], x[i, ])))
rownames(genes.cors) <- rownames(ori)



df.cors_mse <- rbind(genes.mse %>%
  as_tibble() %>%
  mutate(Symbol = rownames(ori)) %>%
  gather(Model, value, 1:4) %>%
  mutate(measure = "mse"),
  genes.cors %>%
    as_tibble() %>%
    mutate(Symbol = rownames(ori)) %>%
    gather(Model, value, 1:4)%>%
    mutate(measure = "cor")
) %>%
  spread(measure, value) %>%
  mutate(Model = factor(Model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense"))) %>%
  left_join(kegg.genes.N, by = "Symbol") %>%
  left_join(kegg.genes.cats, by = "Symbol") %>%
  # left_join(select(annot, gene_type, Symbol), by = "Symbol") %>%
  as_tibble() %>%
  mutate(N = ifelse(is.na(N), 0, N),
          path = ifelse(N == 0, "out", "in"))



png("figures/genes_mse_vs_cor_models.png")
ggplot(df.cors_mse, aes(x = mse, y = cor)) +
  geom_point() +
  facet_wrap(~ Model) +
  theme_bw()
dev.off()

#
# plot(genes.mse, genes.cors)
#
# summary(lm(genes.cors ~exp(genes.mse)))
# points(genes.mse, exp(genes.mse)*-0.4467695 + 1.4435025,  col = "blue") ## Función para pasar de mse a cor
#
# summary(lm(mse ~ path, df.cors_mse ))
# tapply( mse.df$mse,  mse.df$path, summary)
# summary(lm(mse ~ N, mse.df))
#
# tapply( mse.df$mse,  mse.df$Metabolism > 0, summary)
# tapply( mse.df$mse,  mse.df$Cellular > 0, summary)
# tapply( mse.df$mse,  mse.df$Environment > 0, summary)
# tapply( mse.df$mse,  mse.df$Genetic > 0, summary)
# tapply( mse.df$mse,  mse.df$Disease > 0, summary)
# tapply( mse.df$mse,  mse.df$Organism > 0, summary)

png("figures/mse_Npathway.png")
df.cors_mse  %>%
  mutate(class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
         class = factor(class, levels = c("0", "1", "2", "3-5", "6+"))) %>%
  gather(measure, Value, 3:4) %>%
  ggplot(aes(x = class, y = Value)) +
  geom_boxplot() +
  xlab("N pathways per gene") +
  theme_bw() +
  facet_grid(Model ~ measure)
dev.off()


png("figures/mse_category.png", width = 1000)
df.cors_mse  %>%
  filter(N > 0) %>%
  select(-N, -path) %>%
  gather(Category, N, 5:10) %>%
  filter(N > 0) %>%
  gather(measure, Value, 3:4) %>%
  ggplot(aes(x = Category, y = Value)) +
  geom_boxplot() +
  facet_grid(Model ~ measure) +
  theme_bw()
dev.off()

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

inds.mses <- sapply(mat_list, function(x) colMeans((x - ori)**2))


df.inds_mse <- inds.mses %>%
  as_tibble() %>%
  mutate(barcode  = colnames(vst),
         tumor = labels) %>%
  gather(Model, mse, 1:4) %>%
  mutate(Model = factor(Model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense"))) %>%
  left_join(phenos, by = "barcode")


png("figures/mse_individuals.png", width = 3000)
df.inds_mse %>%
  ggplot(aes(x = tumor, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()



png("figures/mse_individuals_windsor.png", width = 3000)
df.inds_mse %>%
  mutate(mse = ifelse(mse > 2, 2, mse)) %>%
  ggplot(aes(x = tumor, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()



png("figures/mse_individuals_condition.png")
df.inds_mse %>%
  mutate(mse = ifelse(mse > 2, 2, mse)) %>%
  ggplot(aes(x = condition, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()

png("figures/mse_individuals_seqcenter.png")
df.inds_mse %>%
  mutate(mse = ifelse(mse > 2, 2, mse)) %>%
  ggplot(aes(x = center, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()

inds.mses.kegg <- sapply(mat_list, function(x) colMeans((x[rownames(x) %in% kegg.map$Symbol,] - ori[rownames(ori) %in% kegg.map$Symbol,])**2))
df.inds_mse_kegg <- inds.mses.kegg %>%
  as_tibble() %>%
  mutate(barcode  = colnames(vst),
         tumor = labels) %>%
  gather(Model, mse, 1:4) %>%
  mutate(Model = factor(Model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense"))) %>%
  left_join(phenos, by = "barcode")


png("figures/mse_individuals_kegg_genes.png", width = 3000)
df.inds_mse_kegg %>%
  ggplot(aes(x = tumor, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()

#
# tapply( inds.mse,  labels, summary) %>% data.frame()
# tapply( inds.mse4,  labels, summary) %>% data.frame()
# tapply( inds.mse5,  labels, summary) %>% data.frame()
#
# ori.ctl <- ori[, labels == "Normal"]
# mat3.ctl <- mat.v3[, labels == "Normal"]
# ori.laml <- ori[, labels == "TCGA-LAML"]
# mat3.laml <- mat.v3[, labels == "TCGA-LAML"]
#
# genes.mse.ctl <- rowMeans((mat3.ctl - ori.ctl)**2)
# genes.mse.laml <- rowMeans((mat3.laml - ori.laml)**2)
#
#
# mats.v2 <- lapply(10:12, function(i){
#   mat <- read_table(paste0("results/TCGA_gexp_norm/kegg_filt_v2.", i, "/model_features/autoencoder_output.tsv.gz"))
#   mat <- mat %>% data.matrix() %>% t()
#   rownames(mat) <- as.character(genes$V1)
#   mat
# })
#
#
# genes.mse.v2 <- lapply(mats.v2, function(x) rowMeans((x - ori)**2))


## From here, old!!
##################################################################################################
## Load kegg - primed
paths <- read.table("results/TCGA_gexp_kegg/v1.1/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])

primed <- readPathways("kegg_v3.1", letters[1:5], paths.vec)

primed.cors <- sapply(paths.vec, pathwayCorr, path_list = primed)
colnames(primed.cors) <- paths.vec
rownames(primed.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")

df.primed <- makeDFsum(primed.cors, "Dropout 0%")

## Get dfs for all dropouts
cors.drops <- lapply(6:9, function(x){
  mod <- paste0("kegg_v2.", x)
  drop20 <- readPathways(mod)

  cors20 <- sapply(paths.vec, pathwayCorr, path_list = drop20)
  colnames(cors20) <- paths.vec
  rownames(cors20 ) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")

  cors20
})


## Get dfs for all dropouts
df.drops <- Map(makeDFsum, cors.drops, paste("Dropout", c("20%", "50%", "75%", "90%")))
df.drops$cot <- df.primed
df.all <- Reduce(rbind, df.drops ) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path")


png("figures/minCor_pathwayCat.png", width = 1000)
df.all %>%
ggplot(aes(x = top_cat, y = minCor, color = model)) +
geom_boxplot() +
theme_bw()
dev.off()

png("figures/minCor_pathwayCat2.png", width = 3000, height = 1000)
df.all %>%
ggplot(aes(x = category, y = minCor, color = model)) +
geom_boxplot() +
theme_bw() +
facet_wrap(~top_cat, scales = "free_x", nrow = 3)
dev.off()


df.all %>% select(-class) %>% spread(model, minCor) %>%
select(-path) %>% data.matrix() %>%
heatmap(scale = "none")

df.all %>% select(-minCor) %>% spread(model, class) %>%
select(-path) %>% data.matrix() %>%
heatmap(scale = "none")

df.sum <- df.all %>% group_by(path) %>% summarize(nlow = sum(class == "low"), nhigh = sum(class == "High"))
goodPaths <- subset(df.sum, nlow == 0)$path


## Filter pathways kegg
paths.filt <- read.table("results/TCGA_gexp_norm/kegg_filt_v3.1/model_trained/pathways_names.txt", header = TRUE)
paths.filt <- as.character(paths.filt[, 1])

primed.filt <- readPathways("kegg_filt_v3.1", letters[1:5], paths.filt)

primed.filt.cors <- sapply(paths.filt, pathwayCorr, path_list = primed.filt)
colnames(primed.filt.cors) <- paths.filt
rownames(primed.filt.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")

df.primed.filt <- makeDFsum(primed.filt.cors, "Dropout 0%")
a <- left_join(df.primed.filt, df.primed, by = "path")
table(filtered = a$class.x, original = a$class.y)

## Get dfs for all dropouts
cors.drops <- lapply(6:9, function(x){
  mod <- paste0("kegg_v2.", x)
  drop20 <- readPathways(mod)

  cors20 <- sapply(paths.vec, pathwayCorr, path_list = drop20)
  colnames(cors20) <- paths.vec
  rownames(cors20 ) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")

  cors20
})


## Get dfs for all dropouts
df.drops <- Map(makeDFsum, cors.drops, paste("Dropout", c("20%", "50%", "75%", "90%")))
df.drops$cot <- df.primed
df.all <- Reduce(rbind, df.drops ) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path")



## hipathia
paths.hip <- read.table("results/TCGA_gexp_norm/hipathia_v3.1/model_trained/pathways_names.txt", header = TRUE, sep = "\t")
paths.hip <- c(as.character(paths.hip[, 1]), "hsa")


## Get dfs for all dropouts
hip.cors.drops <- lapply(6:9, function(x){
  mod <- paste0("hipathia_v2.", x)
  drop20 <- readPathways(mod, letters[1:5], paths.hip)

  cors20 <- sapply(paths.hip, pathwayCorr, path_list = drop20)
  colnames(cors20) <- paths.hip
  rownames(cors20 ) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")

  cors20
})
drop0 <- readPathways("hipathia_v3.1", letters[1:5], paths.hip)
cors0 <- sapply(paths.hip, pathwayCorr, path_list = drop0)
colnames(cors0) <- paths.hip
rownames(cors0 ) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
df.0 <- makeDFsum(cors0, "Dropout 0%")


## Get dfs for all dropouts
df.drops.hip <- Map(makeDFsum, hip.cors.drops, paste("Dropout", c("20%", "50%", "75%", "90%")))

df.drops.hip$cot <- df.0
df.hip.all <- Reduce(rbind, df.drops.hip ) %>%
  mutate(pathway = gsub("-[0-9 ]*$", "", path))

df.sum.hip <- df.hip.all %>% group_by(path, pathway) %>% summarize(nlow = sum(class == "low"), nhigh = sum(class == "High"))
table(low = df.sum.hip$nlow, high =  df.sum.hip$nhigh)


df.sum.hip2 <- df.sum.hip %>% group_by(pathway) %>%
  summarize(nHigh = sum(nhigh > 0), nlow = sum(nlow > 2), nPaths = n())
filter(df.sum.hip2, gsub("P-", "", pathway) %in% cancer.keg)
