##################################################################################
##################################################################################
# CNN network for autoencoder
## Trial 1
##################################################################################
##################################################################################

from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv2D, MaxPooling2D, Dense, Activation, Flatten, Input

def model_generator(x_train, y_train):
    ## Define function to be passed
    def model_function(filters, kernel_size, stride, dense_layer_sizes1, dense_layer_sizes2, pool, alpha):
      input_shape = (x_train.shape[1], x_train.shape[2], 1) 
      num_classes = len(y_train[0])

      model = Sequential()
      ## *********** First layer Conv
      model.add(Conv2D(filters, kernel_size = (2, kernel_size), strides = stride,
      input_shape = input_shape))
      model.add(Activation('relu'))
      model.add(MaxPooling2D((1, pool)))
      model.add(Flatten())
      model.summary()

      model.add(Dense(dense_layer_sizes1, activation='relu'))
      model.add(Dense(dense_layer_sizes2, activation='relu'))
      model.add(Dense(num_classes, activation='softmax'))
      model.output_shape
      opt = Adam(learning_rate = params.alpha)
      model.compile(loss = 'categorical_crossentropy',
                    optimizer = opt,
                    metrics = ['accuracy'], 
                    )
      model.summary()
      return model
    return model_function

def model_generator_train(x_train, y_train, params):
  input_shape = (x_train.shape[1], x_train.shape[2], 1) 
  num_classes = len(y_train[0])

  model = Sequential()
  ## *********** First layer Conv
  model.add(Conv2D(params.filters, kernel_size = (2, params.kernel_size), strides = params.stride,
  input_shape = input_shape))
  model.add(Activation('relu'))
  model.add(MaxPooling2D((1, params.pool)))
  model.add(Flatten())
  
  model.add(Dense(params.dense_layer_sizes1, activation='relu'))
  model.add(Dense(params.dense_layer_sizes2, activation='relu'))
  model.add(Dense(num_classes, activation='softmax'))
  return model
