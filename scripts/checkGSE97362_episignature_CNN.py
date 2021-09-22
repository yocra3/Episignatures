
#'#################################################################################
#'#################################################################################
#'  Test small models in gse97362 using subsets of features
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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, classification_report
from keras.regularizers import l2, l1, l1l2

label_encoder = LabelEncoder()
onehot_encoder = OneHotEncoder(sparse=False)

# Load selected CpGs for both diseases
## Load methylation
f = h5py.File("results/gse97362_episignature_validation/discov_allFeatsassays.h5", 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

f = h5py.File("results/gse97362_episignature_validation/val_allFeatsassays.h5", 'r')
meth_matrix = f['assay001']
x_test = meth_matrix[...]
f.close()

# Load labels and ccnvert them to integers 
with open('results/gse97362_episignature_validation/discov_labels.txt','r') as file:
    label = file.read()
label = label.split('\n')[0:-1]
int1 = label_encoder.fit_transform(label)
int_dim = int1.reshape(len(int1), 1)
y_train = onehot_encoder.fit_transform(int_dim)

with open('results/gse97362_episignature_validation/val_labels.txt','r') as file:
    label = file.read()
label = label.split('\n')[0:-1]
int1 = label_encoder.fit_transform(label)
int_dim = int1.reshape(len(int1), 1)
y_test = onehot_encoder.fit_transform(int_dim)


input_shape = (x_train.shape[1], 1)
num_classes = len(y_train[0])

# Train model ####
model = Sequential()
model.add(Dense(128, activation='relu', input_dim = x_train.shape[1]))))
model.add(Dense(num_classes, activation='softmax'))

opt = Adam(learning_rate = 1e-3)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
history = model.fit(x_train, y_train,
  batch_size = 16,
  epochs = 200,
  verbose = 1, validation_data = (x_test, y_test))


Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)
confusion_matrix(y_class, y_pred)
## Muy buena predicción, son solo alguna muestra que se pierde
# array([[18,  2,  0],
#        [ 0, 85,  0],
#        [ 0,  0,  8]])



# Load all CpGs for both diseases
## Load methylation
f = h5py.File("results/gse97362_episignature_validation/discov_assays.h5", 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

f = h5py.File("results/gse97362_episignature_validation/val_assays.h5", 'r')
meth_matrix = f['assay001']
x_test = meth_matrix[...]
f.close()

# Train model ####
model = Sequential()
model.add(Dense(128, activation='relu', input_dim = x_train.shape[1], kernel_regularizer=l2(0.01)))))
model.add(Dense(num_classes, activation='softmax'))

opt = Adam(learning_rate = 1e-3)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
history = model.fit(x_train, y_train,
  batch_size = 16,
  epochs = 5000,
  verbose = 1, validation_data = (x_test, y_test))


Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)
confusion_matrix(y_class, y_pred)
# array([[12,  7,  1],
#        [ 5, 71,  9],
#        [ 0,  1,  7]])


# Load selected features for both diseases
## Load methylation
f = h5py.File("results/gse97362_episignature_validation/discov_allFlattenassays.h5", 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

f = h5py.File("results/gse97362_episignature_validation/val_allFlattenassays.h5", 'r')
meth_matrix = f['assay001']
x_test = meth_matrix[...]
f.close()

# Train model ####
model = Sequential()
model.add(Dense(128, activation='relu', input_dim = x_train.shape[1], kernel_regularizer=l2(0.01)))
model.add(Dense(num_classes, activation='softmax'))

opt = Adam(learning_rate = 1e-3)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
history = model.fit(x_train, y_train,
  batch_size = 16,
  epochs = 1000,
  verbose = 1, validation_data = (x_test, y_test))


Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)
confusion_matrix(y_class, y_pred)
## Empeora la predicción. Aunque la mayoría de las muestras están bien clasificadas.
# array([[16,  3,  1],
#        [ 0, 84,  1],
#        [ 2,  1,  5]])



# Load all features for both diseases
## Load methylation
f = h5py.File("results/gse97362_episignature_validation/discov_wholeFlattenassays.h5", 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

f = h5py.File("results/gse97362_episignature_validation/val_wholeFlattenassays.h5", 'r')
meth_matrix = f['assay001']
x_test = meth_matrix[...]
f.close()


# Train model ####
model = Sequential()
model.add(Dense(128, activation='relu', input_dim = x_train.shape[1], kernel_regularizer=l2(0.001)))
model.add(Dense(num_classes, activation='softmax'))

opt = Adam(learning_rate = 1e-3)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
history = model.fit(x_train, y_train,
  batch_size = 16,
  epochs = 500,
  verbose = 1, validation_data = (x_test, y_test))


Y_pred = model.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)
confusion_matrix(y_class, y_pred)
# array([[ 9, 11,  0],
#        [ 4, 81,  0],
#        [ 0,  5,  3]])

