##################################################################################
##################################################################################
# Pruning for DNN gexp network 1
##################################################################################
##################################################################################

import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam
import tensorflow_model_optimization as tfmot


alpha = 0.0001
epochs = 10
batch_size = 128

initial_sparsity = 0.7
final_sparsity = 0.995
begin_step = 70
end_step = 700


def apply_pruning_to_dense(layer):
  pruning_params = {
    'pruning_schedule': tfmot.sparsity.keras.PolynomialDecay(initial_sparsity = initial_sparsity,
                                                             final_sparsity = final_sparsity,
                                                             begin_step = begin_step,
                                                             end_step = end_step)
  }
  if layer.name == "dense":
    return tfmot.sparsity.keras.prune_low_magnitude(layer, **pruning_params)
  return layer

def copy_weights(model1, model2):
    w = model1.get_weights()
    wp = model2.get_weights()
    wp[0] = w[0]
    wp[1] = w[1]
    wp[5] = w[2]
    wp[6] = w[3]
    wp[7] = w[4]
    wp[8] = w[5]
    model2.set_weights(wp)
    return model2

def prune_model(model):

    model_for_pruning = tf.keras.models.clone_model(
        model,
        clone_function = apply_pruning_to_dense,
    )
    model_for_pruning = copy_weights(model, model_for_pruning)
    opt = Adam(learning_rate = alpha)
    model_for_pruning.compile(loss='categorical_crossentropy',
      optimizer = opt,
      metrics = ['categorical_accuracy'])
    return model_for_pruning
