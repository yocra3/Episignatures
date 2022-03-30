docker run -v $PWD:$PWD -it yocra3/episignatures_python:1.4 bash
cd /home/SHARED/PROJECTS/Episignatures
python

#'#################################################################################
#'#################################################################################
#' Combine methy from hdf5 data and labels
#'#################################################################################
#'#################################################################################

import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder
import pickle

# Load and reshape data
f = h5py.File('results/TCGA_gexp_combat/vsd_normassays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

with open('results/TCGA_gexp_combat/individuals_labels.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:-1]



## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('results/TCGA_gexp_combat/assay_reshaped.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()


## Normalize genes
means = np.mean(x_train, axis = 0)
ranges = np.max(x_train, axis = 0) - np.min(x_train, axis = 0)
x_train_mod = (x_train - means)/ranges

# Save reshaped training data
f = h5py.File('results/TCGA_gexp/assay_reshaped_norm.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train_mod.shape[0], x_train_mod.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train_mod
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

pickle.dump( [means, ranges], open( "results/TCGA_gexp/norm_values.pb", "wb" ), protocol = 4 )

stds = np.std(x_train, axis = 0)
x_train_std = (x_train - means)/stds

# Save reshaped training data
f = h5py.File('results/TCGA_gexp_combat/assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train_std.shape[0], x_train_std.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train_std
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

pickle.dump( [means, stds], open( "results/TCGA_gexp_combat/standardized_values.pb", "wb" ), protocol = 4 )



# Repeat for proetin coding filtered data
# Load and reshape data
f = h5py.File('results/TCGA_gexp_combat_coding/vsd_normassays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

with open('results/TCGA_gexp_combat_coding/individuals_labels.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:-1]


## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('results/TCGA_gexp_combat_coding/assay_reshaped.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

means = np.mean(x_train, axis = 0)
stds = np.std(x_train, axis = 0)
x_train_std = (x_train - means)/stds

# Save reshaped training data
f = h5py.File('results/TCGA_gexp_combat_coding/assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train_std.shape[0], x_train_std.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train_std
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

# Repeat for GO filtered data
# Load and reshape data
f = h5py.File('results/TCGA_gexp_go/vsd_normassays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

with open('results/TCGA_gexp_go/individuals_labels.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:-1]



## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('results/TCGA_gexp_go/assay_reshaped.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()


## Go and autosomic genes and labels with sex
# Load and reshape data
f = h5py.File('results/TCGA_gexp_go/vsd_norm_autoassays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

with open('results/TCGA_gexp_go/individuals_labels_sex.txt','r') as file:
    project = file.read()
project = project.split('\n')[0:-1]



## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('results/TCGA_gexp_go/assay_reshaped_autosomic_sex.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train.shape[0], x_train.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()
