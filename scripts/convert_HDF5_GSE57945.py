
#'#################################################################################
#'#################################################################################
#' Convert GSE57945 to HDF5 and subset input genes
#'#################################################################################
#'#################################################################################

import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder

# Load and reshape data
f = h5py.File('results/SRP042228/vsd_norm_TCGAgenes_assays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()


# Save reshaped training data
f = h5py.File('results/SRP042228/assay_reshaped.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_input[...] = x_train
f.close()

# Transform data based on TCGA
A = open('results/TCGA_gexp/norm_values.pb', 'rb')
[means, ranges] = pickle.load(A)
A.close()

x_train_mod = (x_train - means)/ranges


f = h5py.File('results/SRP042228/assay_reshaped_norm.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_input[...] = x_train_mod
f.close()


## Repeat for GO genes
# Load and reshape data
f = h5py.File('results/SRP042228/vsd_norm_GOgenes_assays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()


# Save reshaped training data
f = h5py.File('results/SRP042228/assay_reshaped_GOgenes.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_input[...] = x_train
f.close()

## Repeat for GO genes and autosomics
# Load and reshape data
f = h5py.File('results/SRP042228/vsd_norm_autosom_GOgenes_assays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()


# Save reshaped training data
f = h5py.File('results/SRP042228/assay_reshaped_autosom_GOgenes.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_input[...] = x_train
f.close()
