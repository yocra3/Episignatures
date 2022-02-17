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
--name TCGA_gexp_go_sex --params_name v1.6 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v6.py  -resume -profile docker

## TCGA DNN gexp GO sex - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go_sex --params_name v1.6.prun.1 --step prune --model results/TCGA_gexp_go_sex/v1.6/model_trained/TCGA_gexp_go_sex  \
--network_params conf/network_params/dnn_gexp_network1_pruning_1.py -resume -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5 \
--name TCGA_gexp_go_sex --params_name v1.6.prun.2 --step prune --model results/TCGA_gexp_go_sex/v1.6/model_trained/TCGA_gexp_go_sex  \
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

## SRP042228 DNN gexp go sex - v1.6 pruning 1 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_autosom_GOgenes.h5 \
--name SRP042228 --params_name v1.6.prun.1.sex --step features --model results/TCGA_gexp_go_sex/v1.6.prun.1/model_trained/TCGA_gexp_go_sex  \
-resume -profile docker


## SRP042228 DNN gexp go sex - v1.6 pruning 2 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_autosom_GOgenes.h5 \
--name SRP042228 --params_name v1.6.prun.2.sex --step features --model results/TCGA_gexp_go_sex/v1.6.prun.2/model_trained/TCGA_gexp_go_sex  \
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


## Train TCGA gexp pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/kegg_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name kegg_v3.1 --step features \
--model  results/TCGA_gexp_norm/kegg_v3.1/model_trained/TCGA_gexp_norm -resume -profile docker


## DNN gexp - v1.4
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name v1.4 --step train --network conf/network_params/dnn_gexp_network1.py \
--network_params conf/network_params/params_dnn_gexp_network1_v4.py  -resume -profile docker

## TCGA DNN gexp - v1.4 pruning
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name v1.4.prun.2 --step prune --model results/TCGA_gexp_norm/v1.4/model_trained/TCGA_gexp_norm  \
--network_params conf/network_params/dnn_gexp_network1_pruning_2.py -resume -profile docker

## Train TCGA gexp pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1 --step train_biopathways --probes results/TCGA_gexp/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v1.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/hipathia_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.1 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp/assay_reshaped_norm.h5 \
--name TCGA_gexp_norm --params_name hipathia_v3.1 --step features \
--model  results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm -resume -profile docker

## SRP042228 gexp kegg pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_norm.h5 \
--name SRP042228 --params_name kegg_v3.1 --step features --model  results/TCGA_gexp_norm/kegg_v3.1/model_trained/TCGA_gexp_norm   \
-resume -profile docker

## SRP042228 gexp hipathia pathways - v3.1
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_norm.h5 \
--name SRP042228 --params_name hipathia_v3.1 --step features --model  results/TCGA_gexp_norm/hipathia_v3.1/model_trained/TCGA_gexp_norm   \
-resume -profile docker
