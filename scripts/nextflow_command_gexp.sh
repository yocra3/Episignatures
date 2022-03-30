#'#################################################################################
#'#################################################################################
# Nextflow commands for gene expression
#'#################################################################################
#'#################################################################################


## DNN gexp - v1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v1.py  -resume -profile docker

## Extract features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1 --step features --model  results/TCGA_gexp/v1/model_trained/TCGA_gexp/ \
-resume -profile docker

## DNN gexp - v1.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.2 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v2.py  -resume -profile docker

## DNN gexp - v1.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.3 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v3.py  -resume -profile docker

## DNN gexp - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.4 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## DNN gexp - v1.5
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.5 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v5.py  -resume -profile docker

## DNN gexp - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v3.1 --step train --network conf/network_params/dnn_gexp_network3.py \
--network_params conf/network_params/params_dnn_gexp_network3_v1.py  -resume -profile docker


## DNN autoencode gexp - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name autoencod_v1.4 --step autoencoder --network conf/network_params/dnn_gexp_autoencod_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## DNN autoencode gexp - v1.5
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name autoencod_v1.5 --step autoencoder --network conf/network_params/dnn_gexp_autoencod_network1.py \
--network_params conf/network_params/params_dnn_autoencod_gexp_network1_v5.py  -resume -profile docker

## TCGA DNN gexp - v1.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.4 --step features --model results/TCGA_gexp/v1.4/model_trained/TCGA_gexp  \
-resume -profile docker

## TCGA DNN gexp - v1.5 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.5 --step features --model results/TCGA_gexp/v1.5/model_trained/TCGA_gexp  \
-resume -profile docker

## TCGA DNN gexp - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v3.1 --step features --model results/TCGA_gexp/v3.1/model_trained/TCGA_gexp  \
-resume -profile docker

## TCGA autoencode gexp - v1.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name autoencod_v1.4 --step features --model results/TCGA_gexp/autoencod_v1.4/autoencoder_trained/TCGA_gexp_autoencoder  \
-resume -profile docker

## TCGA autoencode gexp - v1.5 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name autoencod_v1.5 --step features --model results/TCGA_gexp/autoencod_v1.5/autoencoder_trained/TCGA_gexp_autoencoder  \
-resume -profile docker


## TCGA DNN gexp - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.4.prun.1 --step prune --model results/TCGA_gexp/v1.4/model_trained/TCGA_gexp  \
--network_params conf/network_params/dnn_gexp_network1_pruning_1.py -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v1.4.prun.2 --step prune --model results/TCGA_gexp/v1.4/model_trained/TCGA_gexp  \
--network_params conf/network_params/dnn_gexp_network1_pruning_2.py -resume -profile docker


## DNN gexp GO - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_go --params_name v1.4 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## TCGA DNN gexp GO - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_go --params_name v1.4.prun.1 --step prune --model results/TCGA_gexp_go/v1.4/model_trained/TCGA_gexp_go  \
--network_params conf/network_params/dnn_gexp_network1_pruning_1.py -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_go --params_name v1.4.prun.2 --step prune --model results/TCGA_gexp_go/v1.4/model_trained/TCGA_gexp_go  \
--network_params conf/network_params/dnn_gexp_network1_pruning_2.py -resume -profile docker


## DNN autoencode gexp GO - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_go --params_name autoencod_v1.4 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## DNN autoencode gexp GO pruning - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_go --params_name autoencod_v1.4.prun.2 --step prune --autoencoder autoencoder --model results/TCGA_gexp_go/autoencod_v1.4/model_trained/TCGA_gexp_go  \
--network_params conf/network_params/dnn_gexp_network1_autoencod_pruning_2.py -resume -profile docker


## DNN autoencode gexp GO autosomic - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go --params_name autoencod__autosom_v1.4 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## DNN autoencode gexp GO autosomic pruning - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go --params_name autoencod_autosom_v1.4.prun.2 --step prune --autoencoder autoencoder --model results/TCGA_gexp_go/autoencod__autosom_v1.4/model_trained/TCGA_gexp_go  \
--network_params conf/network_params/dnn_gexp_network1_autoencod_pruning_2.py -resume -profile docker

## DNN gexp GO sex - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go_sex --params_name v1.7 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v6.py  -resume -profile docker

## TCGA DNN gexp GO sex - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go_sex --params_name v1.7.prun.1 --step prune --model results/TCGA_gexp_go_sex/v1.7/model_trained/TCGA_gexp_go_sex  \
--network_params conf/network_params/dnn_gexp_network1_pruning_1.py -resume -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go_sex --params_name v1.7.prun.2 --step prune --model results/TCGA_gexp_go_sex/v1.7/model_trained/TCGA_gexp_go_sex  \
--network_params conf/network_params/dnn_gexp_network1_pruning_2.py -resume -profile docker

## TCGA autoencode gexp GO - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_kegg --params_name autoencod.kegg --step features --model results/TCGA_gexp_kegg/model_kegg/  \
-resume -profile docker


## SRP042228 DNN gexp - v1.4 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped.h5 \
--name SRP042228 --params_name v1.4 --step features --model results/TCGA_gexp/v1.4/model_trained/TCGA_gexp  \
-resume -profile docker

## SRP042228 DNN gexp - v1.5 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped.h5 \
--name SRP042228 --params_name v1.5 --step features --model results/TCGA_gexp/v1.5/model_trained/TCGA_gexp  \
-resume -profile docker

## SRP042228 DNN gexp go - v1.4 pruning 2 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_GOgenes.h5 \
--name SRP042228 --params_name v1.4.prun.2 --step features --model results/TCGA_gexp_go/v1.4.prun.2/model_trained/TCGA_gexp_go  \
-resume -profile docker

## SRP042228 DNN gexp go sex - v1.7 pruning 1 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_autosom_GOgenes.h5 \
--name SRP042228 --params_name v1.7.prun.1.sex --step features --model results/TCGA_gexp_go_sex/v1.7.prun.1/model_trained/TCGA_gexp_go_sex  \
-resume -profile docker


## SRP042228 DNN gexp go sex - v1.7 pruning 2 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_autosom_GOgenes.h5 \
--name SRP042228 --params_name v1.7.prun.2.sex --step features --model results/TCGA_gexp_go_sex/v1.7.prun.2/model_trained/TCGA_gexp_go_sex  \
-resume -profile docker

## SRP042228 autoencode gexp - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped.h5 \
--name SRP042228 --params_name autoencod_v1.4 --step features --model results/TCGA_gexp/autoencod_v1.4/autoencoder_trained/TCGA_gexp_autoencoder   \
-resume -profile docker

## SRP042228 DNN gexp autoencode go - v1.4 pruning 2 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_GOgenes.h5 \
--name SRP042228 --params_name autoencod_v1.4.prun.2 --step features --model results/TCGA_gexp_go/autoencod_v1.4.prun.2/model_trained/TCGA_gexp_go  \
-resume -profile docker

## SRP042228 DNN gexp autoencode go sex - v1.4 pruning 2 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_autosom_GOgenes.h5 \
--name SRP042228 --params_name autoencod_autosom_v1.4.prun.2 --step features --model results/TCGA_gexp_go/autoencod_autosom_v1.4.prun.2/model_trained/TCGA_gexp_go  \
-resume -profile docker

## SRP042228 autoencode gexp - v1.5
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped.h5 \
--name SRP042228 --params_name autoencod_v1.5 --step features --model results/TCGA_gexp/autoencod_v1.4/autoencoder_trained/TCGA_gexp_autoencoder   \
-resume -profile docker


## DNN gexp - v2.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped.h5 \
--name TCGA_gexp --params_name v2.1 --step train --network conf/network_params/dnn_gexp_network2.py \
--network_params conf/network_params/params_dnn_gexp_network2_v1.py  -resume -profile docker


## Train TCGA gexp pathways - v1.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.1 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network1_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v1.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.2 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network1_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v1.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.3 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network1_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v2.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.1 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network2.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v2.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.2 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network2.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v2.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.3 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network2.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker


## Train TCGA gexp pathways - v2.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.4 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network2.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## Train TCGA gexp pathways - v2.5
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.5 --step train_biopathways --probes results/TCGA_gexp_go/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v5.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network2.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker


## Train TCGA gexp pathways - v1.1 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.1 --step features --model  results/TCGA_gexp_kegg/v1.1/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v1.2 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.2 --step features --model  results/TCGA_gexp_kegg/v1.2/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v1.3 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v1.3 --step features --model  results/TCGA_gexp_kegg/v1.3/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v2.1 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.1 --step features --model  results/TCGA_gexp_kegg/v2.1/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v2.2 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.2 --step features --model  results/TCGA_gexp_kegg/v2.2/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v2.4 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.4 --step features --model  results/TCGA_gexp_kegg/v2.4/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - v2.5 - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped.h5 \
--name TCGA_gexp_kegg --params_name v2.5 --step features --model  results/TCGA_gexp_kegg/v2.5/model_trained/TCGA_gexp_kegg  \
-resume -profile docker

## Train TCGA gexp pathways - primed - features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5  \
--name TCGA_gexp_kegg --params_name primed --step features --model  results/TCGA_gexp_kegg/kegg_primed  \
-resume -profile docker


## Train TCGA gexp autoencoder v2.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name autoencod_v2.1 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network2.py \
--network_params conf/network_params/params_dnn_gexp_autoencod_network2_v1.py  -resume -profile docker

## Train TCGA gexp autoencoder v2.1 -features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name autoencod_v2.1 --step features \
--model  results/TCGA_gexp_norm/autoencod_v2.1/model_trained/TCGA_gexp_norm -resume -profile docker



## Train TCGA gexp autoencoder coding v2.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v2.2 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network2.py \
--network_params conf/network_params/params_dnn_gexp_autoencod_network2_v2.py  -resume -profile docker

## Train TCGA gexp autoencoder coding v2.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v2.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/autoencod_v2.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp autoencoder coding v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v3.1 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network3.py \
--network_params conf/network_params/params_dnn_gexp_autoencod_network3_v1.py  -resume -profile docker

## Train TCGA gexp autoencoder coding v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v3.1 --step features \
--model  results/TCGA_gexp_combat_coding_std/autoencod_v3.1/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp autoencoder coding v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v3.2 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network3.py \
--network_params conf/network_params/params_dnn_gexp_autoencod_network3_v2.py  -resume -profile docker

## Train TCGA gexp autoencoder coding v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v3.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/autoencod_v3.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp autoencoder coding v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v4.1 --step train --autoencoder autoencoder --network conf/network_params/dnn_gexp_autoencod_network4.py \
--network_params conf/network_params/params_dnn_gexp_autoencod_network4_v1.py  -resume -profile docker

## Train TCGA gexp autoencoder coding v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name autoencod_v4.1 --step features \
--model  results/TCGA_gexp_combat_coding_std/autoencod_v4.1/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v1.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v1.3 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network1_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v1.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v1.3 --step features \
--model  results/TCGA_gexp_norm/kegg_v1.3/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.1 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.1/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.2 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.2 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.2/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.4 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.4 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.4/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.7/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.7 - 5 initializations
for i in {a..e}
do
echo kegg_v2.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7${i} --step features \
--model  results/TCGA_gexp_norm/kegg_v2.7${i}/model_trained/TCGA_gexp_norm -profile docker
done



## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.7/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.7 - 5 initializations
for i in {a..e}
do
echo kegg_v2.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7${i} --step features \
--model  results/TCGA_gexp_norm/kegg_v2.7${i}/model_trained/TCGA_gexp_norm  -profile docker
done


## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7b --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.7b --step features \
--model  results/TCGA_gexp_norm/kegg_v2.7b/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.8
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.8 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.8 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.8/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.8 - 5 initializations
for i in {a..e}
do
echo kegg_v2.8${i}

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.8${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -profile docker

## TCGA gexp pathways - v2.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.8${i} --step features \
--model  results/TCGA_gexp_norm/kegg_v2.8${i}/model_trained/TCGA_gexp_norm  -profile docker
done



## Train TCGA gexp pathways - v2.9
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.9 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.9 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.9 --step features \
--model  results/TCGA_gexp_norm/kegg_v2.9/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.9 - 5 initializations
for i in {a..e}
do
echo kegg_v2.9${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.9${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v2.9${i} --step features \
--model  results/TCGA_gexp_norm/kegg_v2.9${i}/model_trained/TCGA_gexp_norm -profile docker
done


## Train TCGA gexp filtered pathways - v2.10
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.10 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.10 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.10 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.10/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.10 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v2.10${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.10${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.10${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.10${i}/model_trained/TCGA_gexp_norm -profile docker
done

## Train TCGA gexp filtered pathways - v2.11
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.11 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.11 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.11 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.11/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.11 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v2.11${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.11${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.11${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.11${i}/model_trained/TCGA_gexp_norm -profile docker
done



## Train TCGA gexp filtered pathways - v2.12
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.12 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v12.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.12 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.12 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.12/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.12 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v2.12${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.12${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v12.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.12${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.12${i}/model_trained/TCGA_gexp_norm -profile docker
done


## Train TCGA gexp filtered pathways - v2.13
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.13 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v13.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.13 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.13 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.13/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.13 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v2.13${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.13${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v13.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v2.13${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v2.13${i}/model_trained/TCGA_gexp_norm -profile docker
done


## Train TCGA gexp pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1 --step features \
--model  results/TCGA_gexp_norm/kegg_v3.1/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v3.1 - 5 initializations
for i in {a..e}
do
echo kegg_v3.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1${i} --step features \
--model  results/TCGA_gexp_norm/kegg_v3.1${i}/model_trained/TCGA_gexp_norm -profile docker
done



## Train TCGA gexp pathways filtered - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.1 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v3.1/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v3.1 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.1${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.1${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v3.1${i}/model_trained/TCGA_gexp_norm -profile docker
done


## Train TCGA gexp pathways filtered - v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v3.2/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v3.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v3.2${i}/model_trained/TCGA_gexp_norm -profile docker
done

## Train TCGA gexp pathways filtered - v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v4.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v4.1 --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v4.1/model_trained/TCGA_gexp_norm -resume -profile docker




## Train TCGA gexp pathways filtered standardized - v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_standardized.h5 \
--name TCGA_gexp_std --params_name kegg_filt_v3.2 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_standardized.h5 \
--name TCGA_gexp_std --params_name kegg_filt_v3.2 --step features \
--model  results/TCGA_gexp_std/kegg_filt_v3.2/model_trained/TCGA_gexp_std -resume -profile docker

## Train TCGA gexp pathways - v3.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_filt_v3.2${i} --step features \
--model  results/TCGA_gexp_norm/kegg_filt_v3.2${i}/model_trained/TCGA_gexp_norm -profile docker
done




## Train TCGA gexp pathways filtered standardized - v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_standardized.h5 \
--name TCGA_gexp_std --params_name kegg_filt_v4.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_standardized.h5 \
--name TCGA_gexp_std --params_name kegg_filt_v4.1 --step features \
--model  results/TCGA_gexp_std/kegg_filt_v4.1/model_trained/TCGA_gexp_std -resume -profile docker




## Train TCGA gexp pathways filtered combat standardized - v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.2 --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.2 --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v3.2/model_trained/TCGA_gexp_combat_std -resume -profile docker

## Train TCGA gexp pathways - v3.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.2${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.2${i} --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v3.2${i}/model_trained/TCGA_gexp_combat_std -profile docker
done



## Train TCGA gexp pathways filtered combat standardized - v3.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.4 --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.4 --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v3.4/model_trained/TCGA_gexp_combat_std -resume -profile docker

## Train TCGA gexp pathways - v3.4 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.4${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.4${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v3.4${i} --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v3.4${i}/model_trained/TCGA_gexp_combat_std -profile docker
done


## Train TCGA gexp pathways filtered combat standardized - v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v4.1 --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v4.1 --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_std -resume -profile docker


## Train TCGA gexp pathways - v4.1 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v4.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v4.1${i} --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name  kegg_filt_v4.1${i} --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v4.1${i}/model_trained/TCGA_gexp_combat_std -profile docker
done

## Train TCGA gexp pathways filtered combat standardized - v5.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v5.1 --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v5.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v5.1 --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_trained/TCGA_gexp_combat_std -resume -profile docker


## Train TCGA gexp pathways - v5.1 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v5.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name kegg_filt_v5.1${i} --step train_biopathways --probes results/TCGA_gexp_combat/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_std --params_name  kegg_filt_v5.1${i} --step features \
--model  results/TCGA_gexp_combat_std/kegg_filt_v5.1${i}/model_trained/TCGA_gexp_combat_std -profile docker
done




## Train TCGA gexp protein coding pathways filtered combat standardized - v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.2 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt_v3.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp pathways - v3.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.2${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.2${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt_v3.2${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.2 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp pathways - v3.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v3.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.2${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.2${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.2${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done

## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.5
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.5 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v5.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.5 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.5 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.5/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp pathways - v3.5 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v3.5${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.5${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v5.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.5${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.5${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done



## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.6
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.6 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.6 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.6${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.6${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done

## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.7 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.7 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.7/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker

## Train TCGA gexp pathways - v3.7 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v3.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.7${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v3.7${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.7${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done


## Train TCGA gexp pathways - v3.3 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v3.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.3${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v3.3${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt_v3.3${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done



## Train TCGA gexp pathways filtered combat standardized - v4.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v4.1 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v4.1 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v4.1 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v4.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt_v4.1${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt_v4.1${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt_v4.1${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done




## Train TCGA gexp pathways filtered combat standardized - v4.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.2 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v4.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v4.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.2${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v4.2${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.2${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done





## Train TCGA gexp pathways filtered combat standardized - v4.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.3 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.3 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v4.3 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v4.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.3${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v4.3${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done

## Train TCGA gexp pathways filtered combat standardized - v4.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.4 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.4 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.4/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v4.4 - 5 initializations
for i in {a..e}
do
echo kegg_filt_v4.4${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v4.4${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v4.4${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.4${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done


## Train TCGA gexp pathways filtered2 combat coding standardized - v5.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.2 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v5.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v5.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v5.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.2${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v5.2${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.2${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done

## Train TCGA gexp pathways filtered2 combat coding standardized - v5.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.3 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v5.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.3 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v5.3 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v5.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.3${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v5.3${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done



## Train TCGA gexp pathways filtered2 combat coding standardized - v5.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.4 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v5.4 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.4 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.4/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v5.4 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v5.4${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v5.4${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v4.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v5.4${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.4${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done





## Train TCGA gexp pathways filtered2 combat coding standardized - v6.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.1 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v6.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.1 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.1/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v6.1 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v6.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.1${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v6.1${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.1${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done

## Train TCGA gexp pathways filtered2 combat coding standardized - v6.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.2 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v6.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.2 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v6.2 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v6.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.2${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v6.2${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done



## Train TCGA gexp pathways filtered2 combat coding standardized - v6.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.3 --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v6.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.3 --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.3/model_trained/TCGA_gexp_combat_coding_std -resume -profile docker


## Train TCGA gexp pathways - v6.3 - 5 initializations
for i in {a..e}
do
echo kegg_filt2_v6.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name kegg_filt2_v6.3${i} --step train_biopathways --probes results/TCGA_gexp_combat_coding/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/kegg_filt_manual_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_gexp_combat_coding_std --params_name  kegg_filt2_v6.3${i} --step features \
--model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.3${i}/model_trained/TCGA_gexp_combat_coding_std -profile docker
done


## SRP042228 gexp kegg fitlered2 pathways - v6.2
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_coding_std_gse.h5 \
--name SRP042228 --params_name kegg_filt2_v6.2 --step features --model  results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_trained/TCGA_gexp_combat_coding_std   \
-resume -profile docker








## DNN gexp - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name v1.4 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## TCGA DNN gexp - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name v1.4.prun.2 --step prune --model results/TCGA_gexp_norm/v1.4/model_trained/TCGA_gexp_norm  \
--network_params conf/network_params/dnn_gexp_network1_pruning_2.py -resume -profile docker


## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7 --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.7/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.7 - 5 initializations
for i in {a..e}
do
echo hipathia_v2.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7${i} --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.7${i}/model_trained/TCGA_gexp_norm -profile docker
done



## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7 --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.7/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.7 - 5 initializations
for i in {a..e}
do
echo hipathia_v2.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7${i} --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.7${i}/model_trained/TCGA_gexp_norm  -profile docker
done


## Train TCGA gexp pathways - v2.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7b --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.7b --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.7b/model_trained/TCGA_gexp_norm -resume -profile docker


## Train TCGA gexp pathways - v2.8
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.8 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.8 --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.8/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.8 - 5 initializations
for i in {a..e}
do
echo hipathia_v2.8${i}

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.8${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -profile docker

## TCGA gexp pathways - v2.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.8${i} --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.8${i}/model_trained/TCGA_gexp_norm  -profile docker
done



## Train TCGA gexp pathways - v2.9
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.9 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker

## TCGA gexp pathways - v2.9 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.9 --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.9/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v2.8 - 5 initializations
for i in {a..e}
do
echo hipathia_v2.9${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.9${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network2_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network1.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v2.9${i} --step features \
--model  results/TCGA_gexp_norm/hipathia_v2.9${i}/model_trained/TCGA_gexp_norm -profile docker
done

## Train TCGA gexp pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1 --step features \
--model  results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm -resume -profile docker

## Train TCGA gexp pathways - v3.1 - 5 initializations
for i in {a..e}
do
echo hipathia_v3.1${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1${i} --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1${i} --step features \
--model  results/TCGA_gexp_norm/hipathia_v3.1${i}/model_trained/TCGA_gexp_norm -profile docker
done





## SRP042228 gexp kegg pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_norm.h5 \
--name SRP042228 --params_name kegg_v3.1 --step features --model  results/TCGA_gexp_norm/kegg_v3.1/model_trained/TCGA_gexp_norm   \
-resume -profile docker

## SRP042228 gexp hipathia pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_norm.h5 \
--name SRP042228 --params_name hipathia_v3.1 --step features --model  results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm   \
-resume -profile docker
