#! /opt/conda/envs/episignatures/bin/python

#'#################################################################################
#'#################################################################################
#'  Train TCGA network
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

from numpy import array, argmax
from keras.models import Sequential, Model
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint

f = h5py.File('./assay_reshaped.h5', 'r')
methy = f['methy'][...]
f.close()


A = open('model.pb', 'rb')
[model, labels] = pickle.load(A)
A.close()

for i in range(3, len(model.layers) - 1):
  new_model = Model(inputs=model.input, outputs=model.layers[i].output)
  Y_pred = new_model.predict(methy)
  df = pd.DataFrame(Y_pred)
  df.to_csv(model.layers[i].name + '.tsv',  sep = "\t", index = False)
