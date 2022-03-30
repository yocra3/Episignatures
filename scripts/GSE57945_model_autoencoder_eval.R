#'#################################################################################
#'#################################################################################
#' Evaluate autoencoder in GSE57945
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(SummarizedExperiment)
library(tidyverse)
library(rjson)
library(rhdf5)
library(HDF5Array)
library(rtracklayer)


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

paths <- read.table("results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
kegg.df.com <- subset(kegg.df, pathID %in% paths.vec)
kegg.genes.N <- kegg.map %>%
  group_by(Symbol) %>%
  summarize(N = n())
#
# annot <- readGFF("/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz")
# annot <- filter(annot, type == "gene") %>%
#   as_tibble() %>%
#  mutate(Symbol = gsub("\\.[0-9]*_[0-9]*", "", gene_id , perl = TRUE))



## Compute mse
ori <- h5read("results/SRP042228/assay_reshaped_coding_std_gse.h5","gexp")
base <- h5read("results/SRP042228/kegg_filt2_v3.6/model_features/autoencoder_output.h5", "auto")
post <- h5read("results/SRP042228/kegg_filt2_v4.3/model_features/autoencoder_output.h5", "auto")
pre <- h5read("results/SRP042228/kegg_filt2_v6.2/model_features/autoencoder_output.h5", "auto")
pre_post <- h5read("results/SRP042228/kegg_filt2_v5.3/model_features/autoencoder_output.h5", "auto")
mat_list <- list(base, post, pre, pre_post)
names(mat_list) <- c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")

vst <- loadHDF5SummarizedExperiment("results/SRP042228/", prefix = "vsd_norm_TCGA_codingGenes_")

genes <- read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")
rownames(ori) <- rownames(base) <- rownames(pre_post) <- rownames(post) <- rownames(pre) <- as.character(genes$V1)


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


mean.mses <- sapply(mat_list, function(x) mean((x - ori)**2))
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

# summary(lm(genes.cors ~exp(genes.mse)))
# points(genes.mse, exp(genes.mse)*-0.4467695 + 1.4435025,  col = "blue") ## Función para pasar de mse a cor

#
# summary(lm(mse ~ path, mse.df))
# tapply( mse.df$mse,  mse.df$path, summary)
# summary(lm(mse ~ N, mse.df))
#
# tapply( mse.df$mse,  mse.df$Metabolism > 0, summary)
# tapply( mse.df$mse,  mse.df$Cellular > 0, summary)
# tapply( mse.df$mse,  mse.df$Environment > 0, summary)
# tapply( mse.df$mse,  mse.df$Genetic > 0, summary)
# tapply( mse.df$mse,  mse.df$Disease > 0, summary)
# tapply( mse.df$mse,  mse.df$Organism > 0, summary)

png("figures/GSE57945_mse_Npathway.png", width = 1200)
df.cors_mse %>%
  mutate(class = ifelse(N > 5, "6+", ifelse(N > 2, "3-5", N)),
         class = factor(class, levels = c("0", "1", "2", "3-5", "6+")),
         mse = ifelse(mse > 1.6, 1.6, mse)) %>%
   gather(measure, Value, 3:4) %>%
  ggplot(aes(x = class, y = Value)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype  = 'dashed') +
  # geom_hline(yintercept = median(genes.mse), linetype  = 'dashed', color = 'blue') +
  theme_bw() +
  xlab("N pathways per gene") +
  facet_grid(measure ~ Model, scales = "free_y" )
dev.off()

png("figures/GSE57945_mse_category.png", width = 1000)
df.cors_mse %>%
  filter(N > 0) %>%
  select(-N, -path) %>%
  mutate(mse = ifelse(mse > 1.6, 1.6, mse)) %>%
  gather(Category, N, 5:10) %>%
  filter(N > 0) %>%
  gather(measure, Value, 3:4) %>%
  ggplot(aes(x = Category, y = Value)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0.3, 0.5, 0.7, 0.9), linetype  = 'dashed') +
  facet_grid(measure ~ Model, scales = "free_y" ) +
  theme_bw()
dev.off()



# png("figures/GSE57945_mse_gene_type.png", width = 1400)
# mse.df %>%
#   mutate(type = ifelse(gene_type %in% c("protein_coding", "processed_pseudogene", "lncRNA", "unprocessed_pseudogene"), gene_type, "Others"),
#           network = ifelse(N > 0, "Network", "Excluded")) %>%
#   gather(measure, val, 2:3) %>%
#   ggplot(aes(x = type, y = val, col = network)) +
#   geom_boxplot() +
#   geom_hline(yintercept = mean.mse, linetype  = 'dashed') +
#   geom_hline(yintercept = median(genes.mse), linetype  = 'dashed', color = 'blue') +
#   theme_bw() +
#   facet_grid(measure ~ set, scales = "free_y" )
# dev.off()
#


inds.mses <- sapply(mat_list, function(x) colMeans((x - ori)**2))


df.inds_mse <- inds.mses %>%
  as_tibble() %>%
  mutate(id  = colnames(vst),
          project = vst$diagnosis2) %>%
  gather(Model, mse, 1:4) %>%
  mutate(Model = factor(Model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")))


png("figures/GSE57945_mse_individuals.png")
df.inds_mse %>%
  ggplot(aes(x = project, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()

png("figures/GSE57945_mse_individuals_windsor.png")
df.inds_mse %>%
  mutate(mse = ifelse(mse > 2, 2, mse)) %>%
  ggplot(aes(x = project, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()


inds.mses.kegg <- sapply(mat_list, function(x) colMeans((x[rownames(x) %in% kegg.map$Symbol,] - ori[rownames(ori) %in% kegg.map$Symbol,])**2))
df.inds_mse_kegg <- inds.mses.kegg %>%
  as_tibble() %>%
  mutate(id  = colnames(vst),
        project = vst$diagnosis2) %>%
  gather(Model, mse, 1:4) %>%
  mutate(Model = factor(Model, levels = c("Pathway", "Pathway + Dense", "Dense + Pathway", "Dense + Pathway + Dense")))

png("figures/GSE57945_mse_individuals_kegg_genes.png")
df.inds_mse_kegg %>%
  ggplot(aes(x = project, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()



png("figures/GSE57945_mse_individuals_kegg_genes_windsor.png")
df.inds_mse_kegg %>%
  mutate(mse = ifelse(mse > 2, 2, mse)) %>%
  ggplot(aes(x = project, y = mse)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(Model ~ .)
dev.off()
