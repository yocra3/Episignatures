#'#################################################################################
#'#################################################################################
# Define genes in KEGG pathways
#'#################################################################################
#'#################################################################################

library(limma)
library(org.Hs.eg.db)
library(tidyverse)
library(rjson)

tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID, column="ENSEMBL", keytype="ENTREZID")
write.table(tab, file = "results/preprocess/kegg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## Check annotation and overlap between pathways
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

paths <- unique(tab$PathwayID)
names(paths) <- paths
kegg.df.com <- subset(kegg.df, pathID %in% paths)



## Discard too general metabolic pathways
met_paths <- subset(kegg.df.com, category == "Global and overview maps")$pathID
met_paths <- met_paths[met_paths != "path:hsa01210"]

paths <- paths[!paths %in% met_paths]

## Remove "path:hsa00524" - pathway of neomicyin sinthesis with high overlap with glucose metabolism
paths <- paths[!paths %in% "path:hsa00524"]


genes_list <- lapply(paths, function(x) subset(tab, PathwayID == x)$Symbol)
tabs.over <- sapply(paths, function(x) sapply(paths, function(y) mean(genes_list[[x]] %in% genes_list[[y]])))
tabs.over <- t(tabs.over)

## Remove redundant pathways: large pathways found in other subpathways
red.paths <-  names(which(colSums(tabs.over > 0.9) > 1))
tab.filt <- subset(tab, PathwayID %in% paths & !PathwayID %in% red.paths)
write.table(tab.filt, file = "results/preprocess/kegg_filt_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")



### New proposal
## Remove pathways from Human Diseases and too complex pathways (after examination in kegg website)
bad_paths <- c("path:hsa04080", "path:hsa04210", "path:hsa04550", "path:hsa04810",
  "path:hsa04640", "path:hsa04613", "path:hsa04062", "path:04922", "path:hsa04916", "path:hsa04928", "path:hsa04935",
  "path:hsa04270", "path:hsa04722", "path:hsa04360")
bad_cats <- c("Cell motility", "Signal transduction", "Signaling molecules and interaction","Transport and catabolism")
out_paths <- c(bad_paths, subset(kegg.df.com, category %in% bad_cats | top_cat == "Human Diseases")$pathID)

tab.filt2 <- subset(tab.filt, !PathwayID %in% out_paths)
write.table(tab.filt2, file = "results/preprocess/kegg_filt_manual_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
