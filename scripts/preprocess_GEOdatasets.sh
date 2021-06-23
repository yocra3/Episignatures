#'#################################################################################
#'#################################################################################
# Preprocess GEO dataset
#'#################################################################################
#'#################################################################################

## Run in data/GSE90496
bgzip -cd GSE90496_beta.txt.gz | awk -F'\t' '{print $1}' |  bgzip -c > probes.txt.gz
bgzip -cd GSE90496_beta.txt.gz | awk -F'\t' '{for (i = 2; i <= NF; i += 2) printf ("%s%c", $i, i + 2 <= NF ? "\t" : "\n");}' |  bgzip -c > beta_vals.txt.gz
bgzip -cd GSE90496_beta.txt.gz | awk -F'\t' '{for (i = 3; i <= NF; i += 2) printf ("%s%c", $i, i + 2 <= NF ? "\t" : "\n");}' |  bgzip -c > detection_pvals.txt.gz

## Run in data/GSE55763
bgzip -cd GSE55763_normalized_betas.txt.gz | awk -F'\t' '{print $1}' |  bgzip -c > probes.txt.gz
bgzip -cd GSE55763_normalized_betas.txt.gz | awk -F'\t' '{for (i = 2; i <= NF; i += 2) printf ("%s%c", $i, i + 2 <= NF ? "\t" : "\n");}' |  bgzip -c > beta_vals.txt.gz
bgzip -cd beta_vals.txt.gz > betas_header.txt
bgzip -cd GSE55763_normalized_betas.txt.gz | awk -F'\t' '{for (i = 3; i <= NF; i += 2) printf ("%s%c", $i, i + 2 <= NF ? "\t" : "\n");}' |  bgzip -c > detection_pvals.txt.gz
