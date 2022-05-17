#'#################################################################################
#'#################################################################################
# Nextflow commands for gene expression
#'#################################################################################
#'#################################################################################


# Use all samples for training
#'#################################################################################

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


# Exclude prostate cancer in training and use combined pathways
#'#################################################################################


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.6
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.6 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo comb_paths_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.6${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.6${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths_v3.6${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.6
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths2_v3.6 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_filt_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths2_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths2_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo comb_paths2_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths2_v3.6${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths2_v3.6${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths2_v3.6${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done



## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.6
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.6${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done

## TCGA gexp PRAD pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5  \
--name TCGA_gexp_coding_PRAD --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## TCGA gexp PRAD tumor pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_tumor_assay_reshaped_standardized.h5  \
--name TCGA_gexp_coding_PRAD_tumor --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker


## TCGA gexp array pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GSE169038/prad_array_reshaped_standardized.h5  \
--name GSE169038 --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker


## GTEX prostate pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/prostate_reshaped_standardized.h5 \
--name GTEx_prostate --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## GTEX testis pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/testis_reshaped_standardized.h5 \
--name GTEx_testis --params_name comb_paths3_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.7
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.7 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.7 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.7 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.7/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.7 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v3.7${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.7${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.7${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.7${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done



## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.8
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.8 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.8 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths_v3.8/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.8 - 5 initializations
for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.8${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths_v3.8${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths_v3.8${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done

## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.8
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.8 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.8 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.8/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.8 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.8${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.8${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.8${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done



## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.9
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.9 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.9 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.9 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.9/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.9 - 5 initializations
for i in {c..e}
do
echo comb_paths3_v3.9${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.9${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v3.9${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.9${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v4.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v4.3 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v4.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v4.3 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v4.3/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v4.3 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v4.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v4.3${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network4.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v4.3${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v4.3${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v5.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v5.3 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v5.3 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v5.3 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v5.3/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v5.3 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v5.3${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v5.3${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v5.3${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v5.3${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done




## Train TCGA gexp protein coding pathways filtered2 combat standardized - v6.2
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v6.2 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v6.2 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v6.2 --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v6.2/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v6.2 - 5 initializations
for i in {a..e}
do
echo comb_paths3_v6.2${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v6.2${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/preprocess/go_kegg_final_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name comb_paths3_v6.2${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/comb_paths3_v6.2${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v2.3
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name autoencod_v2.3 --step train --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_autoencod_network2_v3.py \
--network conf/network_params/dnn_gexp_autoencod_network2.py   -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name autoencod_v2.3  --step features \
--model  results/TCGA_gexp_coding_noPRAD/autoencod_v2.3/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## PRAD
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_PRAD --params_name autoencod_v2.3  --step features \
--model  results/TCGA_gexp_coding_noPRAD/autoencod_v2.3/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker



# Exclude prostate cancer in training and explore mini GO set
#'#################################################################################


## Train TCGA gexp protein coding pathways filtered2 combat standardized - v3.6
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name mini_paths_v3.6 --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/mini_go_gene_map.tsv  -resume -profile docker


## TCGA gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name mini_paths_v3.6 --step features \
--model  results/TCGA_gexp_coding_noPRAD/mini_paths_v3.6/model_trained/TCGA_gexp_coding_noPRAD -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo mini_paths_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name mini_paths_v3.6${i} --step train_biopathways --probes results/TCGA_gexp_coding_noPRAD/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/mini_go_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5 \
--name TCGA_gexp_coding_noPRAD --params_name mini_paths_v3.6${i}  --step features \
--model  results/TCGA_gexp_coding_noPRAD/mini_paths_v3.6${i}/model_trained/TCGA_gexp_coding_noPRAD -profile docker
done
