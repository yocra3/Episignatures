#'#################################################################################
#'#################################################################################
# Nextflow commands for gene expression
#'#################################################################################
#'#################################################################################


# Train main model in GTEx
## Train with all GO terms -> full training
## Train GTEx gexp protein coding pathways  standardized - v3.11
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker

## GTEx gexp pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11 --step features \
--model  results/GTEx_coding/paths_all_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## Train TCGA gexp pathways - v3.6 - 5 initializations
for i in {a..e}
do
echo paths_all_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_all_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done

## Train with subset of pathways 1
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.6 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.6 --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.6/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.6${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.6${i}  --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.6${i}/model_trained/GTEx_coding -profile docker
done



## Train with subset of pathways 1 - higher epochs of step 2
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.11/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt1_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt1_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_filt1_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done

## Train with subset of pathways 2 - higher epochs of step 2
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt2_full_v3.11${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_v3.11${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11${i}/model_trained/GTEx_coding -profile docker
done



## TCGA gexp all - v3.11 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_coding_all --params_name paths_filt2_full_v3.11 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker


## TCGA all with TCGA model
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all_TCGA --params_name  paths_filt2_full_v3.11 --step features \
--model  results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/TCGA_coding_all -resume -profile docker

for i in {a..e}
do
echo TCGA_v3.6${i}
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_coding_all --params_name paths_filt2_full_v3.11${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v6.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all_TCGA --params_name paths_filt2_full_v3.11${i} --step features \
--model  results/TCGA_coding_all/paths_filt2_full_v3.11${i}/model_trained/TCGA_coding_all -resume -profile docker
done

## TCGA all with TCGA model
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all_retrain --params_name  paths_filt2_full_v3.11 --step features \
--model  results/TCGA_coding_all/paths_filt2_full_v3.11/model_trained/GTEX_transfer  -resume -profile docker


## TCGA gexp all pretraining - v3.8 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_coding_all --params_name paths_filt2_pre_v3.8 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all_TCGA --params_name paths_filt2_pre_v3.8 --step features \
--model  results/TCGA_coding_all/paths_filt2_pre_v3.8/model_trained/TCGA_coding_all -resume -profile docker


## TCGA all with GTEx model
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all --params_name  paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## TCGA all with GTEx model
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5 \
--name TCGA_all --params_name paths_filt2_pre_v3.8 --step features \
--model  results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/GTEx_coding -resume -profile docker

## TCGA gexp PRAD pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5  \
--name GTEx_coding_PRAD --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## TCGA gexp PRAD tumor pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_tumor_assay_reshaped_standardized.h5  \
--name GTEx_coding_PRAD_tumor --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## TCGA gexp PRAD tumor pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/TCGA_gexp_coding_noPRAD/prad_ctrl_assay_reshaped_standardized.h5  \
--name GTEx_coding_PRAD_ctrl --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker


## TCGA gexp array pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GSE169038/prad_array_reshaped_standardized.h5  \
--name GSE169038 --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## SRP042228  pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/SRP042228/assay_reshaped_coding_std_gse.h5  \
--name SRP042228 --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## Ibon  pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file nasir_reshaped_standardized.h5  \
--name Ibon --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker

## GSE54460  pathways - v3.6 features
nextflow run workflows/train_model.nf --hdf5_file results/GSE54460/gse54460_reshaped_standardized.h5  \
--name GSE54460 --params_name paths_filt2_full_v3.11 --step features \
--model  results/GTEx_coding/paths_filt2_full_v3.11/model_trained/GTEx_coding -resume -profile docker


# Train GTEx all paths pretraining
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8 --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.8/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.8${i}  --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.8${i}/model_trained/GTEx_coding -profile docker
done

# Train GTEx selected paths pretraining
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8 --step features \
--model  results/GTEx_coding/paths_filt2_pre_v3.8/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo paths_filt2_pre_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v8.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_pre_v3.8${i}  --step features \
--model  results/GTEx_coding/paths_filt2_pre_v3.8${i}/model_trained/GTEx_coding -profile docker
done



# Train GTEx all paths pretraining
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10 --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.10/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.10${i}
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/preprocess/go_kegg_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_all_pretrain_v3.10${i}  --step features \
--model  results/GTEx_coding/paths_all_pretrain_v3.10${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths pretraining
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10 --step features \
--model  results/GTEx_coding/paths_filt2_unfrozen_v3.10/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo unfrozen ${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v10.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_unfrozen_v3.10${i}  --step features \
--model  results/GTEx_coding/paths_filt2_unfrozen_v3.10${i}/model_trained/GTEx_coding -profile docker
done



# Train GTEx selected paths dropout no primed
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7 --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v7.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_noprime_v3.7${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_noprime_v3.7${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths dropout primed
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9 --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_prime_v3.9/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v9.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_drop_prime_v3.9${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_drop_prime_v3.9${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths dropout primed
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3 --step features \
--model  results/GTEx_coding/paths_filt2_full_postdense_v4.3/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v4.3${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network4_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_postdense_v4.3${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_postdense_v4.3${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths pre-post DNN
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3 --step features \
--model  results/GTEx_coding/paths_filt2_full_prepostdense_v5.3/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo prepostdense_v5.3${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network5_v3.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network5.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv  -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_prepostdense_v5.3${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_prepostdense_v5.3${i}/model_trained/GTEx_coding -profile docker
done


# Train GTEx selected paths pre DNN
nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2 --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2 --step features \
--model  results/GTEx_coding/paths_filt2_full_predense_v6.2/model_trained/GTEx_coding -resume -profile docker

for i in {a..e}
do
echo comb_paths_v3.8${i}
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2${i} --step train_biopathways --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network6_v2.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network6.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker


nextflow run workflows/train_model.nf --hdf5_file results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name paths_filt2_full_predense_v6.2${i}  --step features \
--model  results/GTEx_coding/paths_filt2_full_predense_v6.2${i}/model_trained/GTEx_coding -profile docker
done

## Train GTEx autoencoder without pathways
nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name autoencod_v2.3 --step train --probes results/GTEx/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_autoencod_network2_v3.py \
--network conf/network_params/dnn_gexp_autoencod_network2.py   -resume -profile docker

nextflow run workflows/train_model.nf --hdf5_file  results/GTEx/all_reshaped_standardized.h5 \
--name GTEx_coding --params_name autoencod_v2.3  --step features \
--model  results/GTEx_coding/autoencod_v2.3/model_trained/GTEx_coding -resume -profile docker


## Modelo para Kike
nextflow run workflows/train_model.nf --hdf5_file results/Kike/all_reshaped_standardized.h5 \
--name Kike --params_name main_model --step train_biopathways --probes results/Kike/input_genes.txt \
--autoencoder autoencoder --network_params conf/network_params/params_dnn_gexp_pathway_network3_v11.py \
--network conf/network_params/dnn_gexp_autoencod_pathway_network3.py --cpgmap results/GTEx_coding/go_kegg_filt2_gene_map.tsv -resume -profile docker
