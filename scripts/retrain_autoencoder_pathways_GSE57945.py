docker run -it -v /home/SHARED/PROJECTS/Episignatures:/home/SHARED/PROJECTS/Episignatures -w "$PWD" yocra3/episignatures_python:1.4  /bin/bash

python

#'#################################################################################
#'#################################################################################
#'  Retrain autoencoder models from GSE57945
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
import tensorflow_model_optimization as tfmot
from sklearn.model_selection import train_test_split


## GSE57945
f = h5py.File('results/SRP042228/assay_reshaped_std.h5', 'r')
gexp = f['gexp'][...]
f.close()

x_train, x_test = train_test_split(gexp, test_size = 0.25, random_state = 42)

f = h5py.File('results/SRP042228/vsd_norm_TCGAgenes_assays.h5', 'r')
meth_matrix = f['assay001']
vst = meth_matrix[...]
f.close()

A = open('results/TCGA_gexp_combat/standardized_values.pb', 'rb')
[means, stds] = pickle.load(A)
A.close()

vst2 = (vst*0.978904 -4.362919)
vst2_std = (vst2 - means)/stds

vst_std = (vst - np.mean(vst, axis = 0))/np.std(vst, axis = 0)


f = h5py.File('figures/rnaseq_comp/ctrl_combat_assays.h5', 'r')
meth_matrix = f['assay001']
ctrl = meth_matrix[...]
f.close()

f = h5py.File('figures/rnaseq_comp/gse_combat_assays.h5', 'r')
meth_matrix = f['assay001']
geo = meth_matrix[...]
f.close()

geo_std = ((geo*0.9911953 -3.5775030) - means)/stds
ctrl_std = ((ctrl*0.998703 -0.224783) - means)/stds


f = h5py.File('figures/rnaseq_comp/ctrl_combat2_assays.h5', 'r')
meth_matrix = f['assay001']
ctrl2 = meth_matrix[...]
f.close()

f = h5py.File('figures/rnaseq_comp/gse_combat2_assays.h5', 'r')
meth_matrix = f['assay001']
geo2 = meth_matrix[...]
f.close()

geo_std2 = ((geo2*1.0146808 -1.5781602) - means)/stds
ctrl_std2 = ((ctrl2*1.0407685 -1.6213013) - means)/stds


## V3.2
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v3.2/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(x_test)
np.mean((Y_pred - x_test)**2)
# 346.32928

Y_pred2 = model.predict(vst2_std)
np.mean((Y_pred2 - vst2_std)**2)
x_train, x_test = train_test_split(vst2_std, test_size = 0.25, random_state = 42)

np.mean((model.predict(ctrl_std) - ctrl_std)**2)
np.mean((model.predict(geo_std) - geo_std)**2)

np.mean((model.predict(ctrl_std2) - ctrl_std2)**2)
np.mean((model.predict(geo_std2) - geo_std2)**2)


model.layers[0].trainable = False

callbacks = [tfmot.sparsity.keras.UpdatePruningStep()]


model.fit(x_train, x_train,
  batch_size = 64,
  epochs = 50,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test, x_test), workers = 10)
np.mean((model.predict(x_test) - x_test)**2)
# 1.1440233


np.mean((vst_std - model.predict(vst_std))**2)

Y_pred = model.predict(vst_std)

f = h5py.File('results/SRP042228/kegg_filt_v3.2/model_features/autoencoder_output_std.h5', 'w')
dataset_input = f.create_dataset('ori', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = vst_std

dataset_input = f.create_dataset('auto', (Y_pred.shape[0], Y_pred.shape[1]))
dataset_input[...] = Y_pred
f.close()

# 0.8799316521816258
x_train_std, x_test_std = train_test_split(vst_std, test_size = 0.25, random_state = 42)



model.fit(x_train_std, x_train_std,
  batch_size = 64,
  epochs = 20,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test_std, x_test_std), workers = 10)
np.mean((model.predict(x_test_std) - x_test_std)**2)
# 0.7726581485037965

## V4.1
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(x_test)
np.mean((Y_pred - x_test)**2)
# 339.84958
model.layers[0].trainable = False


model.fit(x_train, x_train,
  batch_size = 64,
  epochs = 300,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test, x_test), workers = 10)
np.mean((model.predict(x_test) - x_test)**2)
# 1.1440233

model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v4.1/model_trained/TCGA_gexp_combat_std/')
model.layers[0].trainable = False

np.mean((vst_std - model.predict(vst_std))**2)
# 0.8923994045899205
model.fit(x_train_std, x_train_std,
  batch_size = 64,
  epochs = 20,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test_std, x_test_std), workers = 10)
np.mean((model.predict(x_test_std) - x_test_std)**2)
# 0.788332620673745

## V5.1
model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_trained/TCGA_gexp_combat_std/')
Y_pred = model.predict(x_test)
np.mean((Y_pred - x_test)**2)
# 382.86496
model.layers[0].trainable = False
model.layers[1].trainable = False


model.fit(x_train, x_train,
  batch_size = 64,
  epochs = 300,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test, x_test), workers = 10)
np.mean((model.predict(x_test) - x_test)**2)
# 1.1449416


model = tf.keras.models.load_model('results/TCGA_gexp_combat_std/kegg_filt_v5.1/model_trained/TCGA_gexp_combat_std/')
model.layers[0].trainable = False
model.layers[1].trainable = False

np.mean((vst_std - model.predict(vst_std))**2)
# 0.9239636132502992
model.fit(x_train_std, x_train_std,
  batch_size = 64,
  epochs = 10,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test_std, x_test_std), workers = 10)
np.mean((model.predict(x_test_std) - x_test_std)**2)
# 0.7759203953662229
