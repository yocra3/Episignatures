#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train a neural network defined in Keras
#'#################################################################################
#'#################################################################################

import pickle
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.regularizers import l2, l1

sys.path.append('./')
import model
import params

name = sys.argv[1]
jobs = int(sys.argv[2])

f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

# Train model ####
model = model.model_generator_train(x_train, y_train, params)
opt = Adam(learning_rate = params.alpha)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 5, verbose = 1)]


history = model.fit(x_train, y_train,
  batch_size = params.batch_size,
  epochs = params.epochs,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test), workers = jobs)

history_dict = history.history
pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
model.save( name )
