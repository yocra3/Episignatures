#'#################################################################################
#'#################################################################################
#' Create pathways from hipathia
#'#################################################################################
#'#################################################################################


library(hipathia)
library(tidyverse)
library(org.Hs.eg.db)

pathways <- load_pathways(species = "hsa")



whole_genes <- lapply(pathways$pathigraphs, function(x){
  ig <- x$graph
  genes_list <- V(ig)$genesList
  names(genes_list) <- V(ig)$name
  genes_list
})
names(whole_genes) <- NULL
whole_genes <- unlist(whole_genes, recursive = FALSE)
whole_genes_df <- data.frame(node = rep(names(whole_genes), lengths(whole_genes)),
            GeneID = unlist(whole_genes)) %>%
            filter(!is.na(GeneID) & GeneID != "/") %>%
            as_tibble()

all_effector <- lapply(pathways$pathigraphs, function(x){
  eff <- lapply(x$effector.subgraphs, function(x){
    names(V(x))
  })
  data.frame(PathwayID = rep(names(x$effector.subgraphs), lengths(eff)),
              node = unlist(eff))
})
all_effector <- Reduce(rbind, all_effector)

hipathia_map <- right_join(all_effector, whole_genes_df, by = "node") %>%
  as_tibble()

hipathia_map$Symbol <- mapIds(org.Hs.eg.db, hipathia_map$GeneID, column="ENSEMBL", keytype="ENTREZID")

tab <- hipathia_map[, c("GeneID", "PathwayID", "Symbol")] %>%
  distinct()
write.table(tab, file = "results/preprocess/hipathia_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
