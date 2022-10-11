docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_rsession:1.3  /bin/bash
R

#'#################################################################################
#'#################################################################################
#' Compare stability of scores
#'#################################################################################
#'#################################################################################

##########################################################################

## Load libraries
library(limma)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(GSVA)
library(HDF5Array)
library(hipathia)
library(org.Hs.eg.db)
library(parallel)
library(rjson)
library(ggVennDiagram)


load("data/tcga_gexp_combat.Rdata")
genes <- read.table("./results/GTEx_coding/input_genes.txt")

path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
path_N <- group_by(path.map, PathwayID) %>% summarize(N = n()) %>% mutate(category = PathwayID)

prad.feat.all <- read.table("results/GTEx_coding_PRAD/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)
prad.feat.ctrl <- read.table("results/GTEx_coding_PRAD_ctrl/paths_filt2_full_v3.11/model_features/prune_low_magnitude_dense.tsv", header = TRUE)

paths <- read.table("results/GTEx_coding/paths_filt2_full_v3.11/model_trained/pathways_names.txt", header = TRUE)
paths.vec <- as.character(paths[, 1])
colnames(prad.feat.ctrl) <- colnames(prad.feat.all) <- paths.vec

## Subset data
prad <- gexp_tcga_combat[genes$V1, gexp_tcga_combat$project_id == "TCGA-PRAD"]
prad.vst <- loadHDF5SummarizedExperiment("results/TCGA_gexp_coding_noPRAD/", prefix = "vsd_norm_prad")

## GO + KEGG
path.cor <- sapply(seq_len(ncol(prad.feat.ctrl)), function(i) cor(prad.feat.all[prad$sample_type == "Solid Tissue Normal", i], prad.feat.ctrl[, i] ))

## GSVA
path_genes <- mclapply(paths.vec, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- paths.vec
gsva.all <- gsva(prad, path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
gsva.ctrl <- gsva(prad[, prad$sample_type == "Solid Tissue Normal"], path_genes, min.sz=5, max.sz=500, kcdf = "Poisson")
gsva.all.ctrl <- gsva.all[, gsva.all$sample_type == "Solid Tissue Normal"]

save(gsva.all, gsva.ctrl, file = "results/TCGA_PRAD/GSVA_allPRAD_values.Rdata")
gsva.cor <- sapply(seq_len(nrow(gsva.all.ctrl)), function(i) cor(t(assay(gsva.all.ctrl[i, ])), t(assay(gsva.ctrl[i,])) ))


## hipathia
hip_pathways <- load_pathways(species = "hsa")

trans_data.all <- translate_data(prad.vst, "hsa")
exp_data.all <- normalize_data(trans_data.all)
hip.all <- hipathia(exp_data.all, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.all_vals <- get_paths_data(hip.all )

trans_data.ctrl <- translate_data(prad.vst[, prad.vst$sample_type == "Solid Tissue Normal" ], "hsa")
exp_data.ctrl <- normalize_data(trans_data.ctrl)
hip.ctrl <- hipathia(exp_data.ctrl, hip_pathways, decompose = FALSE, verbose = TRUE)
hip.ctrl_vals <- get_paths_data(hip.ctrl )

save(hip.all_vals, hip.ctrl_vals, file = "results/TCGA_PRAD/hipathia_allPRAD_values.Rdata")

hip.all.ctrl <- hip.all_vals[, hip.all_vals$sample_type == "Solid Tissue Normal"]
hip.cor <- sapply(seq_len(nrow(hip.all.ctrl)), function(i) cor(t(assay(hip.all.ctrl[i, ])), t(assay(hip.ctrl_vals[i,])) ))

df.cor <- data.frame(cor = c(path.cor, gsva.cor, hip.cor),
  Method = rep(c("NetActivity", "GSVA", "Hipathia"), lengths(list(path.cor, gsva.cor, hip.cor)))) %>%
  mutate(Method = factor(Method, levels = c("GSVA", "Hipathia", "NetActivity")))

#
plot_stab <- ggplot(df.cor, aes(x = Method, y = cor)) +
  geom_boxplot() +
  ylab("Scores correlation") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 20))


png("figures/PRAD_scores_stability.png", height = 240)
plot_stab
dev.off()

summary(lm(cor ~ Method, df.cor))

tapply(df.cor$cor, df.cor$Method, summary)
df.cor %>%
  group_by(Method) %>%
  summarize(m = mean(cor > 0.95, na.rm = TRUE))

## Make panel
load("results/GSE169038/de_genes_results.Rdata")
load("results/GSE169038/GSVA_results.Rdata")
load("results/GSE169038/pathways_results.Rdata")
load("results/GSE169038/hipathia.res.Rdata")



# Comparison between PRAD and GEO
## Pathways
load("results/TCGA_PRAD/pathways_results.Rdata")
comb_paths <- left_join(tab.path_prad, tab.paths_geo, by = "category", suffix = c(".TCGA", ".GEO")) %>%
  as_tibble() %>%
  mutate(Signif = ifelse(adj.P.Val.TCGA < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))

path.plot <- ggplot(comb_paths, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("NetActivity") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb_paths), cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))



png("figures/TCGAvsGEO_logFC.png")
path.plot
dev.off()

## Add GO names
term <- read.table("data/GO_terms/term.txt", sep = "\t", quote = "", comment.char = "", as.is = TRUE)
term <- term[, c(2, 4)]
colnames(term) <- c("Name", "category")



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
  mutate(category = paste0("path:hsa", substring(path, 0, 5)),
          Name = gsub("^[0-9]*  ", "", path))

path_map <- rbind(term, dplyr::select(kegg.df, Name, category))
comb_paths_annot <- left_join(comb_paths, path_map, by = "category") %>%
  dplyr::select(category, Name, ends_with("TCGA"), ends_with("GEO"), Signif) %>%
  dplyr::select(-starts_with("AveExpr"), -starts_with("t"), -starts_with("B"))
write.table(comb_paths_annot, file = "results/GSE169038/paths_results_comb.txt",
  sep = "\t", quote = FALSE, row.names = FALSE )

#

summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths ))
summary(lm(logFC.TCGA ~ logFC.GEO, comb_paths, subset = Signif != "None" ))

cor(comb_paths$logFC.TCGA, comb_paths$logFC.GEO)
# [1] 0.5047529
cor(comb_paths[comb_paths$Signif != "None", ]$logFC.TCGA, comb_paths[comb_paths$Signif != "None", ]$logFC.GEO)
# [1] 0.7034897
cor(comb_paths[comb_paths$Signif == "Both", ]$logFC.TCGA, comb_paths[comb_paths$Signif == "Both", ]$logFC.GEO)
# [1] 0.9121616

table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif), margin = 2)
prop.table(table(sign(comb_paths$logFC.TCGA) == sign(comb_paths$logFC.GEO), comb_paths$Signif != "None"), margin = 2)

png("figures/GSE169038_TCGA_propDE_vs_logFCPath.png", height = 300)
rbind(mutate(tab.path_prad, Dataset = "TCGA"),
    mutate(tab.paths_geo, Dataset = "GEO")) %>%
    mutate(Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
    ggplot(aes(x = DE_prop , y = abs(logFC))) +
    geom_point() +
    scale_x_continuous(name = "Proportion of genes DE") +
    scale_y_continuous(name = "logFC Gene Set (absolute value)") +
    facet_wrap(~ Dataset) +
    theme_bw()
dev.off()



## Compare Gene DE
load("results/TCGA_PRAD/genes_results.Rdata")
res_prad$gene <- rownames(res_prad)
comb.genes <- left_join(data.frame(res_prad), tab.genes_geo , by = "gene")   %>%
  mutate(Signif = ifelse(!is.na(padj) & padj  < 0.05, ifelse(adj.P.Val < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val < 0.05, "GEO", "None"))) %>%
  filter(!is.na(pvalue ) & !is.na(P.Value  ))

gene.plot <- ggplot(comb.genes, aes(x = log2FoldChange, y = logFC, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Genes") +
  xlab("log2FC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.genes), cor(comb.genes$log2FoldChange, comb.genes$logFC)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))

png("figures/TCGAvsGEO_genes_logFC.png")
gene.plot
dev.off()


cor(comb.genes$log2FoldChange, comb.genes$logFC)
# [1] 0.2488165
cor(comb.genes[comb.genes$Signif != "None", ]$log2FoldChange, comb.genes[comb.genes$Signif != "None", ]$logFC, use = "complete")
# [1] 0.3250043

table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif), margin = 2)
prop.table(table(sign(comb.genes$log2FoldChange) == sign(comb.genes$logFC), comb.genes$Signif != "None"), margin = 2)


## GSEA
load("results/TCGA_PRAD/GSVA_results.Rdata")
comb.gsva <- inner_join(tab.gsva_prad, tab.gsva_geo, by = "category", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(adj.P.Val.TCGA) & adj.P.Val.TCGA  < 0.05, ifelse(adj.P.Val.GEO < 0.05, "Both", "TCGA"),
                              ifelse(adj.P.Val.GEO < 0.05, "GEO", "None")))


#
gsva.plot <- ggplot(comb.gsva, aes(x = logFC.TCGA, y = logFC.GEO, col = Signif)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("GSVA") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.gsva), cor(comb.gsva$logFC.TCGA, comb.gsva$logFC.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))

png("figures/TCGAvsGEO_GSVA_logFC.png")
gsva.plot
dev.off()

cor(comb.gsva$logFC.TCGA, comb.gsva$logFC.GEO)
# [1] 0.6390369
cor(comb.gsva[comb.gsva$Signif != "None", ]$logFC.TCGA, comb.gsva[comb.gsva$Signif != "None", ]$logFC.GEO, use = "complete")
# [1] 0.809432
cor(comb.gsva[comb.gsva$Signif == "Both", ]$logFC.TCGA, comb.gsva[comb.gsva$Signif == "Both", ]$logFC.GEO, use = "complete")
# [1] 0.9415681


table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif)
prop.table(table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif), margin = 2)
prop.table(table(sign(comb.gsva$logFC.TCGA) == sign(comb.gsva$logFC.GEO), comb.gsva$Signif != "None"), margin = 2)

## hipathia
load("results/TCGA_PRAD/hipathia.res.Rdata")
comb.hipathia <- inner_join(hip.comp_prad, hip.comp_geo, by = "name", suffix = c(".TCGA", ".GEO"))   %>%
  mutate(Signif = ifelse(!is.na(FDRp.value.TCGA) & FDRp.value.TCGA  < 0.05, ifelse(FDRp.value.GEO < 0.05, "Both", "TCGA"),
                              ifelse(FDRp.value.GEO < 0.05, "GEO", "None")))


#
hip.plot <- ggplot(comb.hipathia, aes(x = statistic.TCGA, y = statistic.GEO, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(name = "Significance", values = c("#004D40", "#1E88E5", "#9E9E9E", "#FFC107")) +
  ggtitle("Hipathia") +
  xlab("U statistic in TCGA-PRAD") +
  ylab("U statistic in GEO-PRAD") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(comb.hipathia), cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 20))

png("figures/TCGAvsGEO_hipathia_stat.png")
hip.plot
dev.off()

cor(comb.hipathia$statistic.TCGA, comb.hipathia$statistic.GEO)
# [1] 0.234694
cor(comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.TCGA, comb.hipathia[comb.hipathia$Signif != "None", ]$statistic.GEO, use = "complete")
# [1] 0.4226619

table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif), margin = 2)
prop.table(table(sign(comb.hipathia$statistic.TCGA) == sign(comb.hipathia$statistic.GEO), comb.hipathia$Signif != "None"), margin = 2)

legend <- get_legend(
  # create some space to the left of the legend
  hip.plot + theme(legend.box.margin = margin(0, 0, 0, 12),
                    text = element_text(size = 25))+
                    guides(color = guide_legend(override.aes = list(size = 8)))
)


png("figures/TCGAvsGEO_panel.png", width = 950, height = 1025)
plot_grid(plot_stab,
  plot_grid(
    plot_grid(gene.plot, gsva.plot, hip.plot + theme(legend.position = "none"), path.plot, ncol = 2, labels = LETTERS[2:5], label_size = 20),
    legend, ncol = 2, rel_widths = c(5, 1)
  )
, ncol = 1, labels = c("A", ""), rel_heights = c(1, 3))
dev.off()


## Compare GSVA vs GO+kegg
# paths.vec2 <- gsub(":", "_",paths.vec)
# cor_measures <- sapply(paths.vec2, function(i) cor( geo.feat.filt[, i], geo_gsva[i, ]))
# names(cor_measures) <- paths.vec

gsva_path_df <- left_join(dplyr::select(comb_paths, category, Signif, starts_with("log")),
                          dplyr::select(comb.gsva, category, Signif, starts_with("log")), by = "category", suffix = c(".path", ".GSVA")) %>%
                          mutate(Signif = ifelse(Signif.path == "Both",
                                                    ifelse(Signif.GSVA == "Both", "Both",
                                                            ifelse(Signif.GSVA == "None", "Both paths - None GSVA", "Both paths - 1 GSVA")),
                                                    ifelse(Signif.path == "None",
                                                            ifelse(Signif.GSVA == "Both", "None paths - Both GSVA",
                                                                    ifelse(Signif.GSVA == "None", "None", "None paths - 1 GSVA")),
                                                              ifelse(Signif.GSVA == "Both", "1 paths - Both GSVA",
                                                                    ifelse(Signif.GSVA == "None", "1 paths - None GSVA", "1 paths - 1 GSVA")))),
                                  Signif.GSVA = factor(Signif.GSVA, levels = c("Both", "GEO", "TCGA", "None")),
                                  Signif.path = factor(Signif.path, levels = c("Both", "GEO", "TCGA", "None"))) %>%
                  left_join(path_N, by = "category")


gsva_path_df$GEO_prop <- sapply(gsva_path_df$category, function(path) {
  sel <- subset(comb.genes, gene %in% subset(path.map, PathwayID == path)$Symbol)
  p <- mean(sign(sel$logFC) == 1)
})
gsva_path_df$TCGA_prop <- sapply(gsva_path_df$category, function(path) {
  sel <- subset(comb.genes, gene %in% subset(path.map, PathwayID == path)$Symbol)
  p <- mean(sign(sel$log2FoldChange) == 1)
})


png("figures/TCGAvsGEO_GSVA_GOmodel_comp.png", width = 1000)
gsva_path_df %>% gather(Dataset, logFC, c(3:4, 6:7)) %>%
  mutate(Method = ifelse(grepl("GSVA", Dataset), "GSVA", "GO + KEGG model"),
        Dataset = ifelse(grepl("TCGA", Dataset), "TCGA", "GEO"),
        Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
        spread(Method, logFC) %>%
    filter(Signif != "None") %>%
    ggplot(aes(x = abs(`GO + KEGG model`), y = abs(GSVA), color = Signif.GSVA, shape = Signif.path)) +
    scale_shape_manual(name = "GO + KEGG model", values = c(19, 7, 3, 17)) +
    scale_color_manual(name = "GSVA", values = c("#004D40", "#1E88E5",  "#FFC107", "#9E9E9E")) +
    geom_point() +
    xlab("absolute logFC in GO + KEGG model") +
    ylab("absolute logFC in GSVA") +
    facet_wrap(~ Dataset, scales = "free") +
    theme_bw()
dev.off()

png("figures/TCGAvsGEO_GSVA_GOmodel_comp_prop.png", width = 1000)
gsva_path_df %>%
    gather(Dataset, Proportion, 11:12) %>%
    mutate(GSVA.Sig = ifelse(Dataset == "GEO_prop",
                        ifelse(Signif.GSVA %in% c("Both", "GEO"), "Significant", "Non-significant"),
                        ifelse(Signif.GSVA %in% c("Both", "TCGA"), "Significant", "Non-significant")),
                      path.Sig = ifelse(Dataset == "GEO_prop",
                        ifelse(Signif.path %in% c("Both", "GEO"), "Significant", "Non-significant"),
                        ifelse(Signif.path %in% c("Both", "TCGA"), "Significant", "Non-significant"))) %>%
    gather(Method, Significance, 13:14) %>%
    mutate(Method = recode(Method, GSVA.Sig = "GSVA", path.Sig = "GO + KEGG model"),
            Dataset = recode(Dataset, GEO_prop = "GEO", TCGA_prop = "TCGA" ),
            Dataset = factor(Dataset, levels = c("TCGA", "GEO"))) %>%
    ggplot(aes(color = Significance, x = Proportion)) +
    geom_density() +
    facet_grid(Method ~ Dataset) +
    theme_bw() +
    xlab("Propotion of genes with higher expression in gleason") +
    geom_vline(xintercept = 0.5)
dev.off()




png("figures/TCGAvsGEO_GSVA_signGenes.png", width = 1000)
gsva_path_df %>% gather(Dataset, Proportion, c(11:12)) %>%
        mutate(Prop = pmax(Proportion, 1 - Proportion)) %>%
        ggplot(aes(y = Prop, x = Signif.GSVA)) +
        geom_boxplot() +
        theme_bw() +
        facet_wrap(~Dataset)
dev.off()

png("figures/TCGAvsGEO_GOKEGG_signGenes.png", width = 1000)
gsva_path_df %>% gather(Dataset, Proportion, c(11:12)) %>%
        mutate(Prop = pmax(Proportion, 1 - Proportion)) %>%
        ggplot(aes(y = Prop, x = Signif.path)) +
        geom_boxplot() +
        theme_bw() +
        facet_wrap(~Dataset)
dev.off()

venn_geo <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both", "GEO"))$category,
                              `GO + KEGG` = subset(gsva_path_df, Signif.path %in% c("Both", "GEO"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("GEO") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")
#
venn_tcga <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both", "TCGA"))$category,
                              `GO + KEGG` = subset(gsva_path_df, Signif.path %in% c("Both", "TCGA"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("TCGA") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")

#
venn_both <- ggVennDiagram(list(GSVA = subset(gsva_path_df, Signif.GSVA %in% c("Both"))$category,
                              `GO + KEGG` = subset(gsva_path_df, Signif.path %in% c("Both"))$category),
                            set_size = 7, label_size = 7) +
            scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
            ggtitle("Both") +
            theme(plot.title = element_text(hjust = 0.5, size = 25),
                  legend.position = "none")

png("figures/TCGAvsGEO_GOKEGG_overlap.png", width = 1200, height = 300)
plot_grid(venn_tcga, venn_geo, venn_both, nrow = 1, labels = LETTERS[1:3])
dev.off()
