## Adapt HDF5 Kike
import h5py
import numpy as np
from sklearn.preprocessing import LabelEncoder
import pickle

### Adapt gtex
f = h5py.File('results/Kike/vst_all_assays.h5', 'r')
meth_matrix = f['assay001']
gtex = meth_matrix[...]
f.close()


with open('results/Kike/individuals_labels.txt','r') as file:
    tissue = file.read()
tissue = tissue.split('\n')[0:-1]

## Convert labels to integers labels
# integer encode
label_encoder = LabelEncoder()
label_int = label_encoder.fit_transform(tissue)


means_gtex = np.mean(gtex, axis = 0)
stds_gtex = np.std(gtex, axis = 0)
x_gtex = (gtex - means_gtex)/stds_gtex

# Save reshaped training data
f = h5py.File('results/Kike/all_reshaped_standardized.h5', 'w')
dataset_input = f.create_dataset('methy', (x_gtex.shape[0], x_gtex.shape[1]))
dataset_input[...] = x_gtex
dataset_label = f.create_dataset('label', (len(label_int),))
dataset_label[...] = label_int
dataset_label.attrs['labels'] = label_encoder.classes_.tolist()
f.close()
