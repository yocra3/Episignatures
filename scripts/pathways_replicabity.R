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
library(cowplot)



readPathways <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense.tsv")
    tab <- read.table(path, header = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}


readPathways2 <- function(model, sufix, path_name){
  lapply(sufix, function(i){
    path <- paste0("results/GTEx_coding/", model, i, "/model_features/prune_low_magnitude_dense_1.tsv")
    tab <- read.table(path, header = TRUE)
    tab <- data.matrix(tab)
    colnames(tab) <- path_name
    tab
  })
}


pathwayCorr <- function(path_list, col){
  path_mat <- sapply(path_list, function(x) x[, col])
  cors <- cor(path_mat, method = "spearman")
  cors[upper.tri(cors)]
}

makeDFsum <- function(cors, mod_name){
  df <- data.frame(path = colnames(cors), minCor = colMins(abs(cors)), medCor = colMedians(abs(cors))) %>%
    mutate(class = ifelse(minCor > 0.8, "High", ifelse(minCor < 0.3, "low", "intermediate")),
            model = mod_name) %>%
            as_tibble()

}

kegg.map <- read.table("results/preprocess/go_kegg_gene_map.tsv", header = TRUE)
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

paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.ini <- read.table("results/GTEx_coding/paths_all_full_v3.11/model_trained/pathways_names.txt", header = TRUE)

paths.vec <- as.character(paths[, 1])
paths.ini <- as.character(paths.ini[, 1])

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

all_train <- readCors("paths_all_full_v3.11", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1 + step 2 + step 3")
#
all_nofrozen <- readCors("paths_all_pretrain_v3.10", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1 + step 3")

all_init <- readCors("paths_all_pretrain_v3.8", c("", letters[1:5]), paths.ini, "All gene sets") %>%
  mutate(training = "Step 1")



main <- readCors("paths_filt2_full_v3.11", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1 + step 2 + step 3")
drop <- readCors("paths_filt2_full_drop_noprime_v3.7", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dropout")
pre <- readCors("paths_filt2_pre_v3.8", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1")
drop_full <- readCors("paths_filt2_full_drop_prime_v3.9", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 2 + dropout")

unfrozen <- readCors("paths_filt2_unfrozen_v3.10", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Step 1 + step 3")


# post <- readCors("kegg_filt2_v4.2", c("", letters[1:5]), paths.vec, "Pathway + Dense")  %>%
#   mutate(training = "primed")
post2 <- readCors("paths_filt2_full_postdense_v4.3", c("", letters[1:5]), paths.vec, "Selected gene sets")  %>%
  mutate(training = "Gene Set + Dense")
# post3 <- readCors("kegg_filt2_v4.4", c("", letters[1:5]), paths.vec, "Pathway + Dense")  %>%
#   mutate(training = "primed + dropout")


readCors2 <- function(base, models, paths.name, model.name){
  base <- readPathways2(base, models, paths.name)
  base.cors <- sapply(paths.name, pathwayCorr, path_list = base)
  colnames(base.cors) <- paths.name
  # rownames(base.cors) <- c("1-2", "1-3", "2-3", "1-4", "2-4", "3-4", "1-5", "2-5", "3-5", "4-5")
  df.base <- makeDFsum(base.cors, model.name) %>%
    left_join(data.frame(kegg.N) %>% mutate(path = Var1) %>% select(-Var1), by = "path")

}
# pre <- readCors2("kegg_filt2_v6.1", c("", letters[1:5]), paths.vec, "Dense + Pathway") %>%
#   mutate(training = "primed")
pre2 <- readCors2("paths_filt2_full_predense_v6.2",c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dense + Gene Set")
# pre3 <- readCors2("kegg_filt2_v6.3", c("", letters[1:5]), paths.vec, "Dense + Pathway") %>%
#   mutate(training = "primed + dropout")

# pre.post <- readCors2("kegg_filt2_v5.2", c("", letters[1:5]), paths.vec, "Dense + Pathway + Dense") %>%
#   mutate(training = "primed")
pre2.post <- readCors2("paths_filt2_full_prepostdense_v5.3", c("", letters[1:5]), paths.vec, "Selected gene sets") %>%
  mutate(training = "Dense + Gene Set + Dense")
# pre3.post <- readCors2("kegg_filt2_v5.4", c("", letters[1:5]), paths.vec, "Dense + Pathway + Dense") %>%
#   mutate(training = "primed + dropout")

df.path_sel <- Reduce(rbind, list(main, pre, unfrozen, all_train, all_nofrozen, all_init)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(group = ifelse(model == "Selected gene sets", "Selected GOs and KEGGs", "All GOs + KEGGs"),
         training = recode(training, "pretrained" = "Step 1", "pretrained only" = "Step 1",
                            "whole training" = "Step 1 + step 2 + step 3", "primed + pretrained" = "Step 1 + step 2 + step 3",
                            "unfrozen" = "Step 1 + step 3", "Whole training, without adaptation" = "Step 1 + step 3"),
         training = factor(training, levels = c("Step 1", "Step 1 + step 3",  "Step 1 + step 2 + step 3")))

plot_rep <- ggplot(df.path_sel, aes(x = training, y = minCor)) +
 geom_boxplot() +
 scale_x_discrete(name = "") +
 scale_y_continuous(name = "Replicability") +
 theme_bw() +
 facet_wrap(~ group, scales = "free_x")

png("figures/minCor_pretraning_comp.png", height = 300, width = 800)
plot_rep
dev.off()

plot_genes <-  df.path_sel %>%
  filter(group == "All GOs + KEGGs") %>%
  ggplot(aes(x = Freq, y = minCor)) +
    geom_point() +
    scale_x_log10(name = "Genes per pathway") +
    scale_y_continuous(name = "Replicability") +
    theme_bw() +
    facet_wrap(~ training) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "grey")

png("figures/minCor_pretraning_Ngenes.png", height = 300)
plot_genes
dev.off()

png("figures/replicability_panel.png", height = 600, width = 800)
plot_grid(plot_rep, plot_genes, ncol = 1, labels = c("A", "B"))
dev.off()

table(ifelse(grepl("GO", main$path), "GO", "KEGG"))
df.path_sel %>% group_by(training, group) %>%
summarize(p = mean(minCor > 0.7))

## Models comparison
df.mod <- Reduce(rbind, list(main,post2, pre2, pre2.post)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate( training = recode(training, `Step 1 + step 2 + step 3` = "Gene Set"),
          training = factor(training , levels = c("Gene Set", "Gene Set + Dense", "Dense + Gene Set", "Dense + Gene Set + Dense")))



png("figures/minCor_models.png", width = 800, height = 300)
df.mod %>%
  ggplot(aes(x = training, y = minCor)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Replicability") +
  xlab("Network structure")
dev.off()

## Training comparison
df.train <- Reduce(rbind, list(main, drop, drop_full)) %>%
  left_join(mutate(kegg.df.com, path = pathID) %>% select(path, top_cat, category), by = "path") %>%
  mutate(Training = recode(training, "Step 1 + step 2 + step 3" = "Whole training",
        "Step 2 + dropout" = "Whole training + dropout",
        Dropout = "Step 1 + step 3 + dropout"),
        Training = factor(Training , levels = c("Whole training", "Whole training + dropout", "Step 1 + step 3 + dropout")),)



png("figures/minCor_training.png", width = 600, height = 300)
df.train %>%
  ggplot(aes(x = Training, y = minCor)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Replicability")
dev.off()



#
#
#
# png("figures/minCor_pathwayCat.png",width = 1000)
# df.all %>%
#   filter(training ==  "primed + pretrained")   %>%
#   filter(!is.na(top_cat)) %>%
#   ggplot(aes(x = top_cat, y = minCor, color = model)) +
#   geom_boxplot() +
#   theme_bw() +
#   geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
# dev.off()
#
# png("figures/minCor_pathwayCat2.png", width = 3000, height = 1000)
# df.all %>%
#   filter(training ==  "primed + pretrained")   %>%
#   filter(!is.na(category)) %>%
#   filter(top_cat != "Environmental Information Processing") %>%
#   ggplot(aes(x = category, y = minCor, color = model)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_wrap(~top_cat, scales = "free_x") +
#   geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
# dev.off()
#
# png("figures/minCor_vs_pathwayNgenes.png", width = 800, height = 800)
# df.all %>%
#   filter(training ==  "primed + pretrained" & !is.na(top_cat))   %>%
#   ggplot(aes(y = minCor, x = Freq, color = top_cat)) +
#   geom_point() +
#   scale_x_log10() +
#   theme_bw() +
#   facet_wrap(~model) +
#   geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype = "dashed")
# dev.off()
#
#
# png("figures/medCor_models.png", width = 500)
# df.all %>%
#   ggplot(aes(x = model, y = medCor, color = model)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ training)
# dev.off()


## Plot correlation of worse path
path_vals <- readPathways("paths_filt3_full_v3.6", sufix = c("", letters[1:5]), path_name = paths.vec)
path_mat <- sapply(path_vals, function(m) m[, 427 ])
path_mat2 <- path_mat
path_mat2[, 5] <-  - path_mat2[, 5]

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
