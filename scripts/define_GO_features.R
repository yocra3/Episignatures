#'#################################################################################
#'#################################################################################
#' Define GO terms to include in the model
#'#################################################################################
#'#################################################################################

library(GOfuncR)
library(tidyverse)
library(parallel)


## Load new GO graph
godir <-  "data/GO_terms/"
term <- read.table(paste0(godir, "/term.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)
graph <- read.table(paste0(godir, "/graph_path.txt"), sep = "\t", quote = "", comment.char = "", as.is = TRUE)

genes <- as.character(read.table("./results/TCGA_gexp_combat_coding/input_genes.txt")$V1)


## Get all GO terms
all_gos <- get_child_nodes("GO:0008150", term, graph)

## Remove obsolete terms
term_filt <- subset(term, V5 == 0)
all_gos <- subset(all_gos, child_go_id %in% term_filt$V4)
genes_pairs <- get_anno_genes(go_ids = all_gos$child_go_id, term_df = term, graph_path_df = graph)

tab <- left_join(mutate(all_gos, go_id = child_go_id), genes_pairs, by = "go_id") %>%
  as_tibble() %>%
  filter(!is.na(gene))

tab$Symbol <- mapIds(org.Hs.eg.db, tab$gene , column= "ENSEMBL", keytype="SYMBOL")
tab$PathwayID <- tab$go_id
tab_filt <- subset(tab, Symbol %in% genes)


## Discard GOs with too few or too many genes
gos_gene_tab <- tab_filt %>%
  group_by(go_id) %>%
  summarize(n = n())

bad_gos <- subset(gos_gene_tab,  n < 10 )$go_id

# gos_vec <- unique(gos_df_filt$go_id)
# names(gos_vec) <- gos_vec
# genes_list <- lapply(gos_vec, function(x) subset(gos_df_filt, go_id == x)$gene)
#
# map_mat <- sapply(gos_vec, function(x) sapply(gos_vec, function(y){
#   mean(genes_list[[x]] %in% genes_list[[y]])
# }))
#
#
#
# last_go_vec <- mclapply(unique(gos_df_filt$go_id), function(go){
#   nodes <- get_child_nodes(go, term, graph)
#   !any(subset(nodes, distance == 1)$child_go_id %in% gos_df_filt$go_id)
# }, mc.cores = 20)
# last_go_vec <- unlist(last_go_vec)
#
# good_gos <- unique(gos_df_filt$go_id)[last_go_vec]

# tab <- subset(gos_df_filt, go_id %in% good_gos)
tab_final <- filter(tab_filt, !go_id %in% bad_gos)


write.table(tab[, c("PathwayID", "Symbol")], file = "results/preprocess/go_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

tab2 <- read.table(file = "results/preprocess/kegg_filt_manual_gene_map.tsv", header = TRUE)
tab_com <- rbind(tab_final[, c("PathwayID", "Symbol")], tab2[, c("PathwayID", "Symbol")])

write.table(tab_com, file = "results/preprocess/go_kegg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


gos <- c("GO:0001776", "GO:0002260", "GO:0001782", "GO:0043029")
tab <- subset(gos_df_filt, go_id %in% gos)
tab$Symbol <- mapIds(org.Hs.eg.db, tab$gene , column= "ENSEMBL", keytype="SYMBOL")
tab$GeneID <- mapIds(org.Hs.eg.db, tab$gene, column= "ENTREZID", keytype="SYMBOL")
tab$PathwayID <- tab$go_id
write.table(tab[, c("GeneID", "PathwayID", "Symbol")], file = "results/preprocess/mini_go_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


## Use more stringent threshold
# bad_gos2 <- subset(gos_gene_tab, n > 200 | n < 5)$go_id
# gos_df_filt2 <- filter(gos_df, !go_id %in% bad_gos2)


last_go_vec2 <- mclapply(unique(gos_df_filt$go_id), function(go){
  nodes <- get_child_nodes(go, term, graph)
  sum(subset(nodes, distance == 1)$child_go_id %in% gos_df_filt$go_id) < 2
}, mc.cores = 20)
last_go_vec2 <- unlist(last_go_vec2)
good_gos2 <- unique(gos_df_filt$go_id)[last_go_vec2]
gos_df_filt2 <- filter(gos_df_filt, go_id %in% good_gos2)

last_go_vec3 <- mclapply(unique(gos_df_filt2$go_id), function(go){
  nodes <- get_parent_nodes(go, term, graph)
  !any(subset(nodes, distance == 1)$parent_go_id %in% gos_df_filt2$go_id)
}, mc.cores = 20)
last_go_vec3 <- unlist(last_go_vec3)
good_gos3 <- unique(gos_df_filt2$go_id)[last_go_vec3]

tab2 <- subset(gos_df_filt2, go_id %in% good_gos3)
