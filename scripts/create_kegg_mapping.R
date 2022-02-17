#'#################################################################################
#'#################################################################################
# Define genes in KEGG pathways
#'#################################################################################
#'#################################################################################

library(limma)
library(org.Hs.eg.db)
tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID, column="ENSEMBL", keytype="ENTREZID")
write.table(tab, file = "results/preprocess/kegg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
