#'#################################################################################
#'#################################################################################
#' Compare standardization of scores
#'#################################################################################
#'#################################################################################


## Load libraries
library(tidyverse)
library(cowplot)
library(NetActivityData)
library(DESeq2)
library(parallel)
library(HDF5Array)
library(hipathia)


## Load data
load("results/TCGA_PRAD/pathways_results.Rdata")
load("results/TCGA_PRAD/genes_results.Rdata")
load("results/TCGA_PRAD/GSVA_results.Rdata")
load("results/TCGA_PRAD/hipathia.res.Rdata")

load("results/TCGA_BRCA/pathways_results.Rdata")
load("results/TCGA_BRCA/genes_results.Rdata")
load("results/TCGA_BRCA/GSVA_results.Rdata")
load("results/TCGA_BRCA/hipathia.res.Rdata")



## Compare Gene DE
res_prad$gene <- rownames(res_prad)
res_brca$gene <- rownames(res_brca)

combTum.genes <- left_join(data.frame(res_prad), data.frame(res_brca), by = "gene", suffix = c(".prad", ".brca"))   %>%
  mutate(Signif = ifelse(!is.na(padj.prad) & padj.prad  < 0.05, ifelse(!is.na(padj.brca) & padj.brca < 0.05, "Both", "PRAD"),
                              ifelse(!is.na(padj.brca) & padj.brca < 0.05, "BRCA", "None"))) %>%
  filter(!is.na(pvalue.prad) & !is.na(pvalue.brca))

geneTum.plot <- ggplot(combTum.genes, aes(x = log2FoldChange.prad, y = log2FoldChange.brca, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1BCC78", "#9E9E9E", "#FFC107")) +
  ggtitle("Genes") +
  xlab("log2FC in TCGA-PRAD") +
  ylab("log2FC in TCGA-BRCA") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(combTum.genes), cor(combTum.genes$log2FoldChange.prad, combTum.genes$log2FoldChange.brca)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))


### NetActivity
combTum_path <- left_join(tab.path_brca, tab.path_prad, by = "category", suffix = c(".brca", ".prad")) %>%
  as_tibble() %>%
  mutate(Signif = ifelse(adj.P.Val.brca < 0.05, ifelse(adj.P.Val.prad < 0.05, "Both", "BRCA"),
                              ifelse(adj.P.Val.prad < 0.05, "PRAD", "None")))


#
pathTum.plot <- ggplot(combTum_path, aes(x = logFC.prad, y = logFC.brca, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1BCC78", "#9E9E9E", "#FFC107")) +
  ggtitle("NetActivity") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in TCGA-BRCA") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(combTum_path), cor(combTum_path$logFC.prad, combTum_path$logFC.brca)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))

#
path.map <- read.table("results/GTEx_coding/go_kegg_filt2_gene_map.tsv", header = TRUE)
path_genes <- mclapply(combTum_path$category, function(x) subset(path.map, PathwayID == x & !is.na(Symbol))$Symbol, mc.cores = 10)
names(path_genes) <- combTum_path$category


combTum_path_cord <- filter(combTum_path, sign(logFC.brca) == sign(logFC.prad)) %>%   filter(adj.P.Val.brca < 0.05 & adj.P.Val.prad < 0.05 )
combTum_path_discord <- filter(combTum_path, sign(logFC.brca) != sign(logFC.prad))  %>%   filter(adj.P.Val.brca < 0.05 & adj.P.Val.prad < 0.05 )

data(gtex_gokegg)
weight_cord <- lapply(combTum_path_cord$category, function(go){
  w <- gtex_gokegg[go, ]
  w[w != 0]
})
getCor <- function(w, n){
  genes <- names(tail(sort(abs(w)), n))
  br <- res_brca[genes  , , drop = FALSE]
  pr <- res_prad[genes  , , drop = FALSE]
  colnames(pr) <- paste0("PRAD_", colnames(pr))
  combTum <- cbind(pr, br) %>% as_tibble() %>% dplyr::select(ends_with(c("log2FoldChange", "pvalue", "padj")))
  cor(combTum$PRAD_log2FoldChange, combTum$log2FoldChange)
}
getCor(weight_cord[[1]], 3)
getCor(weight_cord[[2]], 8)
getCor(weight_cord[[3]], 7)
getCor(weight_cord[[4]], 5)
weight_discord <- lapply(combTum_path_discord$category, function(go){
  w <- gtex_gokegg[go, ]
  w[w != 0]
})
getCor(weight_discord[[1]], 7)


cot <- lapply(weight, function(w){
  genes <- names(w)
  br <- res_brca[genes  , , drop = FALSE]
  pr <- res_prad[genes  , , drop = FALSE]
  colnames(pr) <- paste0("PRAD_", colnames(pr))
  comb <- cbind(pr, br) %>% as_tibble() %>% dplyr::select(ends_with(c("log2FoldChange", "pvalue", "padj")))
  comb$gene <- genes
  comb
})

## GSVA

### Compare with PRAD
combTum_gsva <- left_join(tab.gsva_brca, tab.gsva_prad, by = "category", suffix = c(".brca", ".prad")) %>%
  as_tibble()  %>%
  mutate(Signif = ifelse(adj.P.Val.brca < 0.05, ifelse(adj.P.Val.prad < 0.05, "Both", "BRCA"),
                              ifelse(adj.P.Val.prad < 0.05, "PRAD", "None")))

#
gsvaTum.plot <- ggplot(combTum_gsva, aes(x = logFC.prad, y = logFC.brca, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(values = c("#004D40", "#1BCC78", "#9E9E9E", "#FFC107")) +
  ggtitle("GSVA") +
  xlab("logFC in TCGA-PRAD") +
  ylab("logFC in TCGA-BRCA") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(combTum_gsva), cor(combTum_gsva$logFC.prad, combTum_gsva$logFC.brca)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
    text = element_text(size = 20))


comb_gsva_cord <- filter(combTum_gsva, sign(logFC.brca) == sign(logFC.prad)) %>%  filter(adj.P.Val.brca < 0.05 & adj.P.Val.prad < 0.05 )

comb_gsva_discord <- filter(combTum_gsva, sign(logFC.brca) != sign(logFC.prad)) %>%  filter(adj.P.Val.brca < 0.05 & adj.P.Val.prad < 0.05 )



GO_gene_cors_cord <- sapply(comb_gsva_cord$category, function(go){
  br <- res_brca[path_genes[[go]]  , ]
  pr <- res_prad[path_genes[[go]]  , ]
  colnames(pr) <- paste0("PRAD_", colnames(pr))
  comb <- cbind(pr, br) %>% as_tibble() %>% dplyr::select(ends_with(c("log2FoldChange", "pvalue", "padj")))
  cor(comb$PRAD_log2FoldChange, comb$log2FoldChange)
})

GO_gene_cors_discord <- sapply(comb_gsva_discord$category, function(go){
  br <- res_brca[path_genes[[go]]  , ]
  pr <- res_prad[path_genes[[go]]  , ]
  colnames(pr) <- paste0("PRAD_", colnames(pr))
  comb <- cbind(pr, br) %>% as_tibble() %>% dplyr::select(ends_with(c("log2FoldChange", "pvalue", "padj")))
  cor(comb$PRAD_log2FoldChange, comb$log2FoldChange)
})

png("figures/PRAD_BRCA_correlation.png", height = 240)
data.frame(cor = c(GO_gene_cors_cord, GO_gene_cors_discord),
          Type = rep(c("Coherent", "Opposite"), c(length(GO_gene_cors_cord), length(GO_gene_cors_discord)))) %>%
    ggplot(aes(x = cor, fill = Type)) +
          geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.05) +
          theme_bw() +
          xlim(c(-1.1, 1.1)) +
          geom_vline(xintercept = c(-0.7,  0.7), linetype = "dotted") +
          geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed") +
          xlab("PRAD-BRCA effect correlation") +
          ylab("Count")
dev.off()




## hipathia


### Compare with PRAD
combTum_hipathia <- left_join(hip.comp_brca, hip.comp_prad, by = "name", suffix = c(".brca", ".prad")) %>%
  as_tibble()  %>%
  mutate(Signif = ifelse(FDRp.value.brca < 0.05, ifelse(FDRp.value.prad < 0.05, "Both", "TCGA-BRCA"),
                              ifelse(FDRp.value.prad < 0.05, "TCGA-PRAD", "None")),
         Signif = factor(Signif, levels = c("Both", "TCGA-BRCA", "TCGA-PRAD", "GEO-PRAD", "None")))

combTum_hipathia %>% filter(FDRp.value.brca  < 0.05 & FDRp.value.prad  < 0.05 )

#
hipTum.plot <- ggplot(combTum_hipathia, aes(x = statistic.prad, y = statistic.brca, col = Signif)) +
  geom_point() +
  theme_bw()  +
  scale_color_manual(name = "Significance", values = c("#004D40", "#1BCC78",  "#FFC107", "#1E88E5",  "#9E9E9E"), drop = FALSE) +
  ggtitle("Hipathia") +
  xlab("U statistic in TCGA-PRAD") +
  ylab("U statistic in TCGA-BRCA") +
  geom_text(data =  data.frame(label = sprintf("N = %d \n r = %.2f", nrow(combTum_hipathia), cor(combTum_hipathia$statistic.prad, combTum_hipathia$statistic.brca)),
   x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5), aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), col = "black", size = 6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 20))



legendComb <- get_legend(
  # create some space to the left of the legend
  hipTum.plot + theme(legend.box.margin = margin(0, 0, 0, 12),
                    text = element_text(size = 25))+
                    guides(color = guide_legend(override.aes = list(size = 8)))
)

title1 <- ggdraw() +
  draw_label("TCGA vs GEO (PRAD)", size = 30,
                   # element = "text",
                   hjust = 0.5, vjust = 0.5)

#
title2 <- ggdraw() +
  draw_label("PRAD vs BRCA (TCGA)", size = 30,
                   # element = "text",
                   hjust = 0.5, vjust = 0.5)

## Load plots from GEO_PRAD_results.R
png("figures/PRAD_BRCA_GEO_panel.png", width = 1000, height = 1500)
plot_grid(
    plot_grid(title1, gene.plot, gsva.plot, hip.plot + theme(legend.position = "none"), path.plot, ncol = 1, labels = c("", LETTERS[1:4]), label_size = 20, rel_heights = c(0.2, rep(1, 4)), scale = 0.95),
  plot_grid(
    plot_grid(title2, geneTum.plot, gsvaTum.plot, hipTum.plot + theme(legend.position = "none"), pathTum.plot, ncol = 1, labels = c("", LETTERS[5:8]), label_size = 20, rel_heights = c(0.2, rep(1, 4)),  scale = 0.95),
    legendComb, ncol = 2, rel_widths = c(2, 1)
  ), ncol = 2, rel_widths = c(2, 3)
)
dev.off()
