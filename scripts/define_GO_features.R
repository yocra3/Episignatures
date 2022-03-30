#'#################################################################################
#'#################################################################################
#' Define GO terms to include in the model
#'#################################################################################
#'#################################################################################

library(GOfuncR)
library(tidyverse)
library(parallel)

all_gos <- get_child_nodes("GO:0008150")
genes_pairs <- get_anno_genes(go_ids = all_gos$child_go_id)

gos_df <- left_join(mutate(all_gos, go_id = child_go_id), genes_pairs, by = "go_id") %>%
  as_tibble()

## Discard GOs with too few or too many genes
gos_gene_tab <- gos_df %>%
  group_by(go_id) %>%
  summarize(n = n())

bad_gos <- subset(gos_gene_tab, n > 200 | n < 5)$go_id
gos_df_filt <- filter(gos_df, !go_id %in% bad_gos)


last_go_vec <- mclapply(unique(gos_df_filt$go_id), function(go){
  nodes <- get_child_nodes(go)
  !any(subset(nodes, distance == 1)$child_go_id %in% gos_df_filt$go_id)
}, mc.cores = 20)
last_go_vec <- unlist(last_go_vec)

good_gos <- unique(gos_df_filt$go_id)[last_go_vec]

tab <- subset(gos_df_filt, go_id %in% good_gos)
tab$Symbol <- mapIds(org.Hs.eg.db, tab$gene , column= "ENSEMBL", keytype="SYMBOL")
tab$GeneID <- mapIds(org.Hs.eg.db, tab$gene, column= "ENTREZID", keytype="SYMBOL")
tab$PathwayID <- tab$go_id

write.table(tab[, c("GeneID", "PathwayID", "Symbol")], file = "results/preprocess/go_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

tab2 <- read.table(file = "results/preprocess/kegg_filt_manual_gene_map.tsv", header = TRUE)
tab_com <- rbind(tab[, c("GeneID", "PathwayID", "Symbol")], tab2)

write.table(tab_com, file = "results/preprocess/go_kegg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
