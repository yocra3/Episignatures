#! /opt/conda/envs/episignatures/bin/python

#'#################################################################################
#'#################################################################################
#'  Train TCGA network
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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint

sys.path.append('./')
import network_config

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

input_shape = (x_train.shape[1], 1)
num_classes = len(y_train[0])

# Train model ####
model = Sequential()
## *********** First layer Conv
model.add(Conv1D(network_config.filters, kernel_size = network_config.kernel_size, 
    strides = network_config.stride, 
    input_shape = input_shape))
model.add(Activation('relu'))
model.add(MaxPooling1D(2))
## ********* Classification layer
model.add(Flatten())
model.add(Dense(network_config.dense_layer_sizes, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))
model.compile(loss='categorical_crossentropy',
  optimizer = 'adam',
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'categorical_accuracy', patience = 3, verbose=0)]
history = model.fit(x_train, y_train,
  batch_size = 128,
  epochs = 20,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test))

history_dict = history.history
pickle.dump( history_dict, open( "history_model.pb", "wb" ), protocol = 4 )
pickle.dump( model, open( "model.pb", "wb" ), protocol = 4 )
