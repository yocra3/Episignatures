docker run -v $PWD:$PWD -it yocra3/episignatures_python:1.4 bash
cd /home/SHARED/PROJECTS/Episignatures
python

#'#################################################################################
#'#################################################################################
# Define autonecoder based on KEGG pathways
#'#################################################################################
#'#################################################################################


import pickle
import h5py
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import tensorflow as tf
import os
import tempfile
from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
import tensorflow_model_optimization as tfmot
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam

## Load Data
A = open('work/8b/d1baab8d18c0d5ee22f63446bafb23/train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('work/8b/d1baab8d18c0d5ee22f63446bafb23/test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()



## Load genes
with open('results/TCGA_gexp_go/input_genes_autosomics.txt','r') as file:
    genes = file.read()
genes = genes.split('\n')[0:-1]

## Select map rows with cpgs in input
with open('results/preprocess/kegg_gene_map.tsv','r') as file:
    gene_map = pd.read_csv(file, delimiter = ' ', index_col = False)


gene_map_filt = gene_map[gene_map['Symbol'].isin(genes)]
pathway_count = gene_map_filt['PathwayID'].value_counts()

gene_dict = dict(zip(genes, range(len(genes))))
gene_map_filt['idx'] = [gene_dict[gene] for gene in gene_map_filt['Symbol']]

## Create list of inputs - pathway
pathway_dict = {}
for g, idx in zip(gene_map_filt['PathwayID'].values, gene_map_filt['idx'].values):
  if g not in pathway_dict:
    pathway_dict[g] = [idx]
  else:
    pathway_dict[g].append(idx)

## Define gene mask ####
ngenes = len(genes)
npaths = len(pathway_dict.keys())

mask_d  = np.zeros((ngenes, npaths), float)

pcol = 0

pathways = pathway_dict.keys()
for path in pathways:
  mask_d[pathway_dict[path], pcol] = 1
  pcol += 1

df = pd.DataFrame(data = list(pathways) )
df.to_csv('results/TCGA_gexp_go_prune/pathways_names.txt',  sep = "\t", index = False)


## Define model
model = Sequential()
model.add(tfmot.sparsity.keras.prune_low_magnitude(Dense(npaths, input_dim = ngenes, activation = 'relu'),
    tfmot.sparsity.keras.ConstantSparsity(0.5, 100000, end_step =  100000, frequency = 100)))
model.add(Dense(ngenes))

## Add gene mask
w = model.get_weights()
w[2] = mask_d
model.set_weights(w)

opt = Adam(learning_rate = 0.001)
callbacks = [tfmot.sparsity.keras.UpdatePruningStep()]


model.compile(loss='mse',
    optimizer = opt,
    metrics = ['mse'])

history = model.fit(x_train, x_train,
  batch_size = 128,
  epochs = 100,
  callbacks = callbacks,
  verbose = 1, validation_data = (x_test, x_test), workers = 10)
model.save( 'results/TCGA_gexp_kegg/model_kegg' )

mask_log = mask_d == 1

opt = Adam(learning_rate = 0.01)

x_train_mod = (x_train - np.mean(x_train, axis = 0))/np.mean(x_train, axis = 0)
x_test_mod = (x_test - np.mean(x_train, axis = 0))/np.mean(x_train, axis = 0)

mini_weights = list()
for i in range(mask_log.shape[1]):
    print(i)
    minimodel = Sequential()
    minimodel.add(Dense(1, input_dim = np.sum(mask_log[:, i])))
    minimodel.add(Dense(np.sum(mask_log[:, i])))
    minimodel.compile(loss='mse',
        optimizer = opt,
        metrics = ['mse'])
    minimodel.fit(x_train_mod[:, mask_log[:, i]], x_train_mod[:, mask_log[:, i]],
        batch_size = 128,
        epochs = 2,
        verbose = 1, validation_data = (x_test_mod[:, mask_log[:, i]], x_test_mod[:, mask_log[:, i]]), workers = 10)
    w = minimodel.get_weights()
    mini_weights.append([w[0], w[1]])
    np.mean((x_test_mod[:, mask_log[:, i]] - np.mean(x_test_mod[:, mask_log[:, i]], axis = 0)) ** 2)

## Convert weights to matrix
new_w  = np.zeros((ngenes, npaths), float)
new_b  = np.zeros(npaths, float)

pcol = 0
for i in range(len(mini_weights)):
  new_w[mask_log[:, i], i] = mini_weights[i][0][:, 0]
  new_b[i] = mini_weights[i][1]
  pcol += 1


## Define model
model = Sequential()
model.add(tfmot.sparsity.keras.prune_low_magnitude(Dense(npaths, input_dim = ngenes, activation = 'relu'),
    tfmot.sparsity.keras.ConstantSparsity(0.5, 100000, end_step =  100000, frequency = 100)))
model.add(Dense(ngenes))

## Add gene mask
w = model.get_weights()
w[0] = new_w
w[1] = new_b
w[2] = mask_d
model.set_weights(w)

new_model = Model(inputs=model.input, outputs=model.layers[0].output)
Y_pred = new_model.predict(x_test_mod)
np.sum(np.mean(Y_pred == 0, axis = 0) < 1)

opt = Adam(learning_rate = 0.001)
callbacks = [tfmot.sparsity.keras.UpdatePruningStep()]


model.compile(loss='mse',
    optimizer = opt,
    metrics = ['mse'])

history = model.fit(x_train_mod, x_train_mod,
  batch_size = 128,
  epochs = 30,
  callbacks = callbacks,
  verbose = 1, validation_data = (x_test_mod, x_test_mod), workers = 10)


model.save('results/TCGA_gexp_kegg/kegg_primed')

## Export features
f = h5py.File('results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5', 'r')
all_mat = f['methy'][...]
f.close()

all_mat_mod = (all_mat - np.mean(x_train, axis = 0))/np.mean(x_train, axis = 0)
new_model = Model(inputs=model.input, outputs=model.layers[0].output)
Y_pred = new_model.predict(all_mat_mod)
np.sum(np.mean(Y_pred == 0, axis = 0) < 1)
df = pd.DataFrame(Y_pred)
df.to_csv("results/TCGA_gexp_kegg/primed/model_features/prune_low_magnitude_dense_2938.tsv",  sep = "\t", index = False)


w0 = np.concatenate([i[0] for i in mini_weights], axis = 1)

minimodel = Sequential()
# minimodel.add(Dense(10, input_dim = np.sum(mask_log[:, i]), activation = 'relu'))
minimodel.add(Dense(1))
minimodel.add(Dense(np.sum(mask_log[:, i])))
minimodel.compile(loss='mse',
    optimizer = opt,
    metrics = ['mse'])
minimodel.fit(x_train_mod[:, mask_log[:, i]], x_train_mod[:, mask_log[:, i]],
    batch_size = 128,
    epochs = 10,
    verbose = 1, validation_data = (x_test_mod[:, mask_log[:, i]], x_test_mod[:, mask_log[:, i]]), workers = 10)
new_model = Model(inputs=model.input, outputs=model.layers[0].output)


np.mean((x_test_mod  - np.mean(x_test_mod , axis = 0)) ** 2)
np.mean((x_test_mod[:, mask_log[:, i]] - np.mean(x_test_mod[:, mask_log[:, i]], axis = 0)) ** 2)
# model = tf.keras.models.load_model('results/TCGA_gexp_go_sex/v1.6.prun.2/model_trained/TCGA_gexp_go_sex')
# w = model.get_weights()
#
# np.savetxt('results/TCGA_gexp_norm/v1.4.prun.2/model_trained/dense_weights.txt.gz', w[2])
