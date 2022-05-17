docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_python:1.4  /bin/bash

python

#'#################################################################################
#'#################################################################################
#'  Extract features from network
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

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model

## Train
f = h5py.File("results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5", 'r')
methy = f['methy'][...]
f.close()

## PRAD
f = h5py.File("results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5", 'r')
prad = f['methy'][...]
f.close()


## GSE57945
f = h5py.File('results/SRP042228/assay_reshaped_coding_std_gse.h5', 'r')
ibd = f['methy'][...]
f.close()


## V3.6
model = tf.keras.models.load_model('results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_trained/TCGA_gexp_coding_noPRAD')
Y_pred = model.predict(methy)

f = h5py.File('results/TCGA_gexp_coding_noPRAD/comb_paths3_v3.6/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()



Y_prad = model.predict(prad)

f = h5py.File('results/TCGA_gexp_coding_PRAD/comb_paths3_v3.6/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_prad.shape[0], Y_prad.shape[1]))
dataset_input[...] = Y_prad
f.close()

Y_ibd = model.predict(ibd)

f = h5py.File('results/SRP042228/comb_paths3_v3.6/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_ibd.shape[0], Y_ibd.shape[1]))
dataset_input[...] = Y_ibd
f.close()

## autoencoder
auto = tf.keras.models.load_model('results/TCGA_gexp_coding_noPRAD/autoencod_v2.3/model_trained/TCGA_gexp_coding_noPRAD')
Y_pred_a = auto.predict(methy)

f = h5py.File('results/TCGA_gexp_coding_noPRAD/autoencod_v2.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred_a.shape[0], Y_pred_a.shape[1]))
dataset_input[...] = Y_pred_a
f.close()



Y_prad_a = auto.predict(prad)

f = h5py.File('results/TCGA_gexp_coding_PRAD/autoencod_v2.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_prad_a.shape[0], Y_prad_a.shape[1]))
dataset_input[...] = Y_prad_a
f.close()


Y_ibd_a = auto.predict(ibd)

f = h5py.File('results/SRP042228/autoencod_v2.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_ibd_a.shape[0], Y_ibd_a.shape[1]))
dataset_input[...] = Y_ibd_a
f.close()















#### Old
f = h5py.File('./results/TCGA_gexp/assay_reshaped_norm.h5', 'r')
methy = f['methy'][...]
f.close()

## V3.2
model = tf.keras.models.load_model('results/TCGA_gexp_norm/kegg_filt_v3.2/model_trained/TCGA_gexp_norm/')
Y_pred = model.predict(methy)
df = pd.DataFrame(Y_pred)
df.to_csv('results/TCGA_gexp_norm/kegg_filt_v3.2/model_features/autoencoder_output.tsv',  sep = "\t", index = False)

for i in range(10, 13):
    model = tf.keras.models.load_model('results/TCGA_gexp_norm/kegg_filt_v2.' + str(i) + '/model_trained/TCGA_gexp_norm/')
    Y_pred = model.predict(methy)
    df = pd.DataFrame(Y_pred)
    df.to_csv('results/TCGA_gexp_norm/kegg_filt_v2.' + str(i) + '/model_features/autoencoder_output.tsv.gz',  sep = "\t", index = False)

## V4.1
model = tf.keras.models.load_model('results/TCGA_gexp_norm/kegg_filt_v4.1/model_trained/TCGA_gexp_norm/')
Y_pred = model.predict(methy)
df = pd.DataFrame(Y_pred)
df.to_csv('results/TCGA_gexp_norm/kegg_filt_v3.2/model_features/autoencoder_output.tsv.gz',  sep = "\t", index = False)

# Standardized
f = h5py.File('./results/TCGA_gexp_combat/assay_reshaped_standardized.h5', 'r')
methy = f['methy'][...]
f.close()


## V3.2
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v3.2/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(methy)

f = h5py.File('results/TCGA_gexp_combat_std/kegg_filt_v3.2/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V4.1
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(methy)

f = h5py.File('results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V5.1
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(methy)

f = h5py.File('results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## GSE57945
f = h5py.File('results/SRP042228/assay_reshaped_std_gse.h5', 'r')
gexp = f['gexp'][...]
f.close()

## V3.2
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v3.2/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(gexp)

os.mkdir('results/SRP042228/kegg_filt_v3.2/')
os.mkdir('results/SRP042228/kegg_filt_v3.2/model_features/')
f = h5py.File('results/SRP042228/kegg_filt_v3.2/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V4.1
os.mkdir('results/SRP042228/kegg_filt_v4.1/')
os.mkdir('results/SRP042228/kegg_filt_v4.1/model_features/')
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(gexp)

f = h5py.File('results/SRP042228/kegg_filt_v4.1/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V5.1
os.mkdir('results/SRP042228/kegg_filt_v5.1/')
os.mkdir('results/SRP042228/kegg_filt_v5.1/model_features/')

model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(gexp)

f = h5py.File('results/SRP042228/kegg_filt_v5.1/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()



## Standardized and protein_coding
f = h5py.File('results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5', 'r')
ori = f['methy'][...]
f.close()


## filt2 - V3.6
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(ori)

f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## filt2 - V4.3
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(ori)

f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## filt2 - V5.3
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(ori)

f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## filt2 - V6.2
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(ori)

f = h5py.File('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## GSE57945
f = h5py.File('results/SRP042228/assay_reshaped_coding_std_gse.h5', 'r')
gexp = f['gexp'][...]
f.close()


## V3.6
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v3.6/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(gexp)

os.mkdir('results/SRP042228/kegg_filt2_v3.6/')
os.mkdir('results/SRP042228/kegg_filt2_v3.6/model_features/')
f = h5py.File('results/SRP042228/kegg_filt2_v3.6/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V4.3
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v4.3/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(gexp)

os.mkdir('results/SRP042228/kegg_filt2_v4.3/')
os.mkdir('results/SRP042228/kegg_filt2_v4.3/model_features/')
f = h5py.File('results/SRP042228/kegg_filt2_v4.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V5.3
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v5.3/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(gexp)

os.mkdir('results/SRP042228/kegg_filt2_v5.3/')
os.mkdir('results/SRP042228/kegg_filt2_v5.3/model_features/')
f = h5py.File('results/SRP042228/kegg_filt2_v5.3/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

## V6.2
model = tf.keras.models.load_model('results/TCGA_gexp_combat_coding_std/kegg_filt2_v6.2/model_trained/TCGA_gexp_combat_coding_std')
Y_pred = model.predict(gexp)

os.mkdir('results/SRP042228/kegg_filt2_v6.2/')
os.mkdir('results/SRP042228/kegg_filt2_v6.2/model_features/')
f = h5py.File('results/SRP042228/kegg_filt2_v6.2/model_features/autoencoder_output.h5', 'w')
dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()
