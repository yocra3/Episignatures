##################################################################################
##################################################################################
# CNN 1D network for autoencoder
## Trial 1
##################################################################################
##################################################################################

from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input

def model_generator(x_train, y_train):
    ## Define function to be passed
    def model_function(filters, kernel_size, stride, dense_layer_sizes1, dense_layer_sizes2, pool, alpha):
      input_shape = (x_train.shape[1], 1) 
      num_classes = len(y_train[0])

      model = Sequential()
      ## *********** First layer Conv
      model.add(Conv1D(filters, kernel_size = kernel_size, strides = stride, input_shape = input_shape))
      model.add(Activation('relu'))
      model.add(MaxPooling1D(pool))
      model.add(Flatten())

      model.add(Dense(dense_layer_sizes1, activation='relu'))
      model.add(Dense(dense_layer_sizes2, activation='relu'))
      model.add(Dense(num_classes, activation='softmax'))
      model.summary()
      opt = Adam(learning_rate = alpha)
      model.compile(loss = 'categorical_crossentropy',
                    optimizer = opt,
                    metrics = ['accuracy'], 
                    )
      model.summary()
      return model
    return model_function

def model_function(filters, kernel_size, stride, dense_layer_sizes1, dense_layer_sizes2, pool, alpha, input_shape, num_classes):

  model = Sequential()
  ## *********** First layer Conv
  model.add(Conv1D(filters, kernel_size = kernel_size, strides = stride, input_shape = input_shape))
  model.add(Activation('relu'))
  model.add(MaxPooling1D(pool))
  model.add(Flatten())
  
  model.add(Dense(dense_layer_sizes1, activation='relu'))
  model.add(Dense(dense_layer_sizes2, activation='relu'))
  model.add(Dense(num_classes, activation='softmax'))
  model.summary()
  opt = Adam(learning_rate = alpha)
  model.compile(loss = 'categorical_crossentropy',
                optimizer = opt,
                metrics = ['accuracy'], 
  )
  return model


def model_generator_train(x_train, y_train, params):
  input_shape = (x_train.shape[1], 1) 
  num_classes = len(y_train[0])

  model = Sequential()
  ## *********** First layer Conv
  model.add(Conv1D(params.filters, kernel_size = params.kernel_size, strides = params.stride,
  input_shape = input_shape))
  model.add(Activation('relu'))
  model.add(MaxPooling1D(params.pool))
  model.add(Flatten())
  
  model.add(Dense(params.dense_layer_sizes1, activation='relu'))
  model.add(Dense(params.dense_layer_sizes2, activation='relu'))
  model.add(Dense(num_classes, activation='softmax'))
  return model
