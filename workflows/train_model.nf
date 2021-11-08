#!/usr/bin/env nextflow
/*
========================================================================================
Train model on TCGA data
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////
date = java.time.LocalDate.now()
params.name = ""
params.test_prop = 0.2
params.step = "train"

params.network_params = "a.txt"
params.network = "a.txt"
params.model = "a.txt"
params.config = "a.txt"
params.labels = "a.txt"
params.probes = "a.txt"

name = params.name
test_prop = params.test_prop
step = params.step
layer = params.layer

ch_input_hdf5 = Channel.fromPath(file("${params.hdf5_file}"))
ch_input_hdf5_se = Channel.fromPath(file("${params.hdf5_se_file}"))

network = file("${params.network}")
network_params = file("${params.network_params}")
model = file("${params.model}")
labels = file("${params.labels}")
probes = file("${params.probes}")


random_config = file("${params.config}")
params_name = "v1"
params_name = "${params.params_name}"


include { PY_DIVIDE_TRAIN_TEST } from '../modules/local/py_divide_training_test/main.nf' addParams( options: [publish_files : ['s':'']])
include { RANDOM_SEARCH } from '../modules/local/random_search/main.nf' addParams( options: [publish_dir: "${name}/${date}/random_search/${params_name}/"])
include { TRAIN_MODEL } from '../modules/local/py_train_model/main.nf' addParams( options: [publish_dir: "${name}/${date}/model_trained/${params_name}/"])
include { PY_EXPORT_RESULTS } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${date}/model_trained/${params_name}/"])
include { EXTRACT_MODEL_FEATURES } from '../modules/local/extract_model_features/main.nf' addParams( options: [publish_dir: "${name}/${date}/model_features/${params_name}/"])
include { TRANSFER_LEARNING } from '../modules/local/transfer_learning/main.nf' addParams( options: [publish_dir: "${name}/${date}/transfer_learning/${params_name}/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_TRANSFER} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${date}/model_trained/${params_name}/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_AUTOENCOD} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${date}/autoencoder_trained/${params_name}/"])
include { PY_EXPORT_RESULTS as PY_EXPORT_RESULTS_AUTOENCOD2} from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "${name}/${date}/autoencoder_trained/${params_name}/"])

include { R_FILTER_MISSINGS_HDF5 } from '../modules/local/r_filter_missings_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_FILTER_INVARIANTCPG_HDF5 } from '../modules/local/r_filter_invariantcpg_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_SUBSTITUE_MISSING_HDF5 } from '../modules/local/r_substitue_missing_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_CPGS_MEDIANS } from '../modules/local/r_cpgs_medians/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_WRITE_PROJECT_LABELS } from '../modules/local/r_write_project_labels/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_SELECT_INPUTCPG_HDF5 } from '../modules/local/r_select_inputcpgs_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])

include { PY_RESHAPE_HDF5_2D } from '../modules/local/py_reshape_hdf5_2D/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { R_CONVERT2D_HDF5 } from '../modules/local/r_convert2D_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { PY_RESHAPE_HDF5 } from '../modules/local/py_reshape_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { PY_COMBINE_METHY_LABEL_HDF5 } from '../modules/local/py_combine_methy_label_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/${name}/${date}"])
include { TRAIN_AUTOENCODER_MODEL } from '../modules/local/py_train_autoencoder_model/main.nf' addParams( options: [publish_dir: "${name}/${date}/autoencoder_trained/${params_name}/"])


workflow  {

  if (step == "preprocessCNN"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_1D(ch_input)
    }
  if (step == "preprocessCNN_trans"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_1D_TRANSF(ch_input, probes)
  }
  if (step == "preprocessDNN"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_DNN(ch_input)
  }
  if (step == "preprocess2D"){
    ch_input = Channel.of("${name}").concat(ch_input_hdf5).concat(ch_input_hdf5_se)
    WF_PREPROCESS_2D(ch_input)
  }

  if (step == "random"){
    WF_RANDOM_SEARCH(ch_input_hdf5, test_prop, network, random_config, name)
  }
  if (step == "train"){
    WF_TRAIN_MODEL(ch_input_hdf5, test_prop, network, network_params, name)
  }
  if (step == "transfer"){
    WF_TRANSFER_MODEL(ch_input_hdf5, test_prop, model, labels, name, layer)
  }
  if (step == "features"){
    EXTRACT_MODEL_FEATURES(ch_input_hdf5, model)
  }
  if (step == "autoencoder"){
    WF_AUTOENCODER_MODEL(ch_input_hdf5, test_prop, model, labels, name)
  }
}


workflow WF_RANDOM_SEARCH {

  take:
  ch_input_hdf5
  test_prop
  network
  random_config
  name

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  RANDOM_SEARCH( PY_DIVIDE_TRAIN_TEST.out.train.collect(), network, random_config, name)

}

workflow WF_TRAIN_MODEL {

  take:
  ch_input_hdf5
  test_prop
  network
  network_params
  name

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  TRAIN_MODEL( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, network, network_params, name )
  PY_EXPORT_RESULTS( ch_input_hdf5, TRAIN_MODEL.out.history, TRAIN_MODEL.out.model, TRAIN_MODEL.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name )

}

workflow WF_TRANSFER_MODEL {

  take:
  ch_input_hdf5
  test_prop
  model
  labels
  name
  layer

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  TRANSFER_LEARNING( ch_input_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, labels, name, layer )
  PY_EXPORT_RESULTS_TRANSFER( ch_input_hdf5, TRANSFER_LEARNING.out.history, TRANSFER_LEARNING.out.model,  TRANSFER_LEARNING.out.labels, PY_DIVIDE_TRAIN_TEST.out.test, name )

}

workflow WF_AUTOENCODER_MODEL {

  take:
  ch_input_hdf5
  test_prop
  model
  labels
  name

  main:
  PY_DIVIDE_TRAIN_TEST ( ch_input_hdf5, test_prop )
  TRAIN_AUTOENCODER_MODEL( PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, name )
  PY_EXPORT_RESULTS_AUTOENCOD( ch_input_hdf5, TRAIN_AUTOENCODER_MODEL.out.old_history, TRAIN_AUTOENCODER_MODEL.out.model_old,  labels, PY_DIVIDE_TRAIN_TEST.out.test, name )

}


workflow WF_PREPROCESS_1D {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_FILTER_INVARIANTCPG_HDF5( R_FILTER_MISSINGS_HDF5.out )
  R_SUBSTITUE_MISSING_HDF5( R_FILTER_INVARIANTCPG_HDF5.out.res )

  R_CPGS_MEDIANS( R_FILTER_INVARIANTCPG_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_RESHAPE_HDF5( R_SUBSTITUE_MISSING_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}

workflow WF_PREPROCESS_DNN {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_FILTER_INVARIANTCPG_HDF5( R_FILTER_MISSINGS_HDF5.out )
  R_SUBSTITUE_MISSING_HDF5( R_FILTER_INVARIANTCPG_HDF5.out.res )

  R_CPGS_MEDIANS( R_FILTER_INVARIANTCPG_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_COMBINE_METHY_LABEL_HDF5( R_SUBSTITUE_MISSING_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}

workflow WF_PREPROCESS_2D {

  take:
  ch_input_hdf5


  main:
  R_FILTER_MISSINGS_HDF5( ch_input_hdf5.collect() )
  R_FILTER_INVARIANTCPG_HDF5( R_FILTER_MISSINGS_HDF5.out )
  R_SUBSTITUE_MISSING_HDF5( R_FILTER_INVARIANTCPG_HDF5.out.res )

  R_CPGS_MEDIANS( R_FILTER_INVARIANTCPG_HDF5.out.res )
  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )
  R_CONVERT2D_HDF5(R_SUBSTITUE_MISSING_HDF5.out.res )
  PY_RESHAPE_HDF5_2D( R_CONVERT2D_HDF5.out.res, R_WRITE_PROJECT_LABELS.out )


}

workflow WF_PREPROCESS_1D_TRANSF {

  take:
  ch_input_hdf5
  probes

  main:
  R_SUBSTITUE_MISSING_HDF5( ch_input_hdf5.collect() )
  R_SELECT_INPUTCPG_HDF5( R_SUBSTITUE_MISSING_HDF5.out.res , probes )

  R_WRITE_PROJECT_LABELS( ch_input_hdf5.collect() )

  PY_RESHAPE_HDF5( R_SELECT_INPUTCPG_HDF5.out.h5, R_WRITE_PROJECT_LABELS.out )

}
