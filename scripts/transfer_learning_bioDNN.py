#'#################################################################################
#'#################################################################################
#'  Add new layer from bioDNN output in blood
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
import tensorflow as tf
import h5py

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report


## Load dense layer data
input_mat = np.genfromtxt('results/GSE97362/2021-11-22/model_features/bioDNN_v3/dense.tsv', 
delimiter='\t', skip_header = 1)
filt_mat = input_mat[:, np.any(input_mat > 0, axis = 0)]

## Load labels
f = h5py.File('results/preprocess/gse97362/2021-07-15/assay_reshaped.h5', 'r')
labels = f['label'][:]
proj_labels = f['label'].attrs['labels']
f.close()

onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = labels.reshape(len(labels), 1)
onehot_labels = onehot_encoder.fit_transform(integer_encoded)

## Split train/test
x_train, x_test, y_train, y_test = train_test_split(filt_mat, onehot_labels,
                                                    stratify = onehot_labels,
                                                    test_size = 0.4, random_state = 42)

## Define model
num_classes = len(y_train[0])

model = Sequential()
model.add(Dense(100, input_dim = x_train.shape[1], activation='relu'))
model.add(Dense(50, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))
opt = Adam(learning_rate = 0.001)
model.compile(loss = 'categorical_crossentropy',
  optimizer = opt,
  metrics = ['accuracy'])
model.summary()

callbacks = [EarlyStopping(monitor = 'val_loss', patience = 10, verbose = 1)]
history = model.fit(x_train, y_train,
  batch_size = 8,
  epochs = 500,
  verbose = 1, validation_data = (x_test, y_test), workers = 5)

Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)

cm = confusion_matrix(y_class, y_pred)
pd.DataFrame(cm, columns = proj_labels)
#    disease state: CHARGE  disease state: Control  disease state: Kabuki
# 0                      7                       7                      2
# 1                      2                      48                      0
# 2                      1                       1                      6
# 

## Load concat layer data
conc_mat = np.genfromtxt('results/GSE97362/2021-11-22/model_features/bioDNN_v3/concatenate.tsv', 
delimiter='\t', skip_header = 1)
filt_conc = conc_mat[:, np.any(conc_mat > 0, axis = 0)]

## Split train/test
x_train2, x_test2, y_train2, y_test2 = train_test_split(filt_conc, onehot_labels,
                                                    stratify = onehot_labels,
                                                    test_size = 0.4, random_state = 42)



model2 = Sequential()
model2.add(Dense(1000, input_dim = x_train2.shape[1], activation='relu'))
model2.add(Dense(100, activation='relu'))
model2.add(Dense(num_classes, activation='softmax'))
opt = Adam(learning_rate = 0.001)
model2.compile(loss = 'categorical_crossentropy',
  optimizer = opt,
  metrics = ['accuracy'])
model2.summary()

callbacks = [EarlyStopping(monitor = 'val_loss', patience = 10, verbose = 1)]
history = model2.fit(x_train2, y_train2,
  batch_size = 8,
  epochs = 100,
  verbose = 1, validation_data = (x_test2, y_test2), workers = 5)

Y_pred2 = model2.predict(x_test2)
y_pred2 = np.argmax(Y_pred2, axis=1)
y_class2 = np.argmax(y_test2, axis=1)

cm2 = confusion_matrix(y_class2, y_pred2)
pd.DataFrame(cm2, columns = proj_labels)
#    disease state: CHARGE  disease state: Control  disease state: Kabuki
# 0                      9                       7                      0
# 1                      0                      50                      0
# 2                      0                       1                      7
