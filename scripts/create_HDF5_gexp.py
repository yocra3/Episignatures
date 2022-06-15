docker run -v $PWD:$PWD -it yocra3/episignatures_python:1.4 bash
cd /home/SHARED/PROJECTS/Episignatures
python

#'#################################################################################
#'#################################################################################
#' Standardize Gene expression data
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


## Standardize genes
stds = np.std(x_train, axis = 0)
means = np.mean(x_train, axis = 0)
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

# Repeat for protein coding filtered data
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



# Use PRAD as validation
# Load and reshape data
f = h5py.File('results/TCGA_gexp_coding_noPRAD/vsd_norm_trainassays.h5', 'r')
meth_matrix = f['assay001']
x_train = meth_matrix[...]
f.close()

with open('results/TCGA_gexp_coding_noPRAD/individuals_labels.txt','r') as file:
    project = file.read()

project = project.split('\n')[0:-1]


## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(project)


# Save reshaped training data
f = h5py.File('results/TCGA_gexp_coding_noPRAD/train_assay_reshaped.h5', 'w')
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
f = h5py.File('results/TCGA_gexp_coding_noPRAD/train_assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_train_std.shape[0], x_train_std.shape[1]))
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_input[...] = x_train_std
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()

pickle.dump( [means, stds], open( "results/TCGA_gexp_coding_noPRAD/standardized_values.pb", "wb" ), protocol = 4 )

### Change PRAD
f = h5py.File('results/TCGA_gexp_coding_noPRAD/vsd_norm_pradassays.h5', 'r')
meth_matrix = f['assay001']
x_prad = meth_matrix[...]
f.close()

means_prad = np.mean(x_prad, axis = 0)
stds_prad = np.std(x_prad, axis = 0)
x_prad_std = (x_prad - means_prad)/stds_prad

# Save reshaped training data
f = h5py.File('results/TCGA_gexp_coding_noPRAD/prad_assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prad_std.shape[0], x_prad_std.shape[1]))
dataset_input[...] = x_prad_std
f.close()


### Tumor
f = h5py.File('results/TCGA_gexp_coding_noPRAD/vsd_norm_prad_tumorassays.h5', 'r')
meth_matrix = f['assay001']
x_prad_tum = meth_matrix[...]
f.close()

means_prad_tum = np.mean(x_prad_tum, axis = 0)
stds_prad_tum = np.std(x_prad_tum, axis = 0)
x_prad_tum_std = (x_prad_tum - means_prad_tum)/stds_prad_tum

# Save reshaped training data
f = h5py.File('results/TCGA_gexp_coding_noPRAD/prad_tumor_assay_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prad_tum_std.shape[0], x_prad_tum_std.shape[1]))
dataset_input[...] = x_prad_tum_std
f.close()

### Change GSE169038
f = h5py.File('results/GSE169038/network_genesassays.h5', 'r')
meth_matrix = f['assay001']
x_array = meth_matrix[...]
f.close()

means_array = np.mean(x_array, axis = 0)
stds_array = np.std(x_array, axis = 0)
x_prad_array = (x_array - means_array)/stds_array
x_prad_array[:, stds_array == 0] = 0

# Save reshaped training data
f = h5py.File('results/GSE169038/prad_array_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prad_array.shape[0], x_prad_array.shape[1]))
dataset_input[...] = x_prad_array
f.close()


### Adapt gtex
f = h5py.File('results/GTEx/vst_all_assays.h5', 'r')
meth_matrix = f['assay001']
gtex = meth_matrix[...]
f.close()

means_gtex = np.mean(gtex, axis = 0)
stds_gtex = np.std(gtex, axis = 0)
x_gtex = (gtex - means_gtex)/stds_gtex

# Save reshaped training data
f = h5py.File('results/GTEx/all_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_gtex.shape[0], x_gtex.shape[1]))
dataset_input[...] = x_gtex
f.close()


### Adapt prostate
f = h5py.File('results/GTEx/vst_prostate_assays.h5', 'r')
meth_matrix = f['assay001']
prostate = meth_matrix[...]
f.close()

means_prost = np.mean(prostate, axis = 0)
stds_prost = np.std(prostate, axis = 0)
x_prost = (prostate - means_prost)/stds_prost

# Save reshaped training data
f = h5py.File('results/GTEx/prostate_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_prost.shape[0], x_prost.shape[1]))
dataset_input[...] = x_prost
f.close()


### Adapt testis
f = h5py.File('results/GTEx/vst_testis_assays.h5', 'r')
meth_matrix = f['assay001']
testis = meth_matrix[...]
f.close()

means_testis = np.mean(testis, axis = 0)
stds_testis = np.std(testis, axis = 0)
x_testis = (testis - means_testis)/stds_testis

# Save reshaped training data
f = h5py.File('results/GTEx/testis_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_testis.shape[0], x_testis.shape[1]))
dataset_input[...] = x_testis
f.close()

### Change HELIX
f = h5py.File('results/HELIX/network_genesassays.h5', 'r')
meth_matrix = f['assay001']
x_array = meth_matrix[...]
f.close()

means_array = np.mean(x_array, axis = 0)
stds_array = np.std(x_array, axis = 0)
x_helix_array = (x_array - means_array)/stds_array
x_helix_array[:, stds_array == 0] = 0

# Save reshaped training data
f = h5py.File('results/HELIX/helix_array_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_helix_array.shape[0], x_helix_array.shape[1]))
dataset_input[...] = x_helix_array
f.close()
