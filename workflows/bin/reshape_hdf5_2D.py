#! /opt/conda/envs/episignatures/bin/python

#'#################################################################################
#'#################################################################################
#' Reshape hdf5 data to add channel dimension
#'#################################################################################
#'#################################################################################

import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder

# Load and reshape data
f = h5py.File('assays.h5', 'r')
meth_matrix = f['methy']
x_train = meth_matrix[...].reshape(meth_matrix.shape[0], meth_matrix.shape[1], meth_matrix.shape[2], 1)
f.close()

with open('labels.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:-1]

## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('./assay_reshaped.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1], x_train.shape[2], 1))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

