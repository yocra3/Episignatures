#'#################################################################################
#'#################################################################################
#' Generate annotation map for first layer of DNN
#'#################################################################################
#'#################################################################################

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(Other)
annot <- Other
annot <- annot[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]

data(Locations)


## Select CpG positions
annot.cpg <- annot[grep("cg", rownames(annot)), ]
allGenes <- unique(unlist(strsplit(annot.cpg$UCSC_RefGene_Name, ";")))
names(allGenes) <- allGenes

annot.cpg$Genes <- strsplit(annot.cpg$UCSC_RefGene_Name, ";")
annot.list.df <- lapply(seq_len(nrow(annot.cpg)), function(i){
  gene <- unique(annot.cpg$Genes[[i]])
  if (length(gene) == 0){
    gene <- "Intergenic"
  }
  data.frame(gene = gene, 
             cpgs = rownames(annot.cpg)[i])
})
annot.df <- data.frame(gene = unlist(sapply(annot.list.df, function(x) x$gene)),
                       cpgs = unlist(sapply(annot.list.df, function(x) x$cpgs)))
genes_tab <- table(annot.df$gene)
small_genes <- names(genes_tab)[genes_tab < 3]
annot.df[annot.df$gene %in% small_genes, "gene"] <- "Intergenic"
annot.df$chromosome <- Locations[annot.df$cpgs, "chr"]

write.table(annot.df, file = "results/preprocess/cpg_gene_map.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

