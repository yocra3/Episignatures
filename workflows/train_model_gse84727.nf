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

hdf5_ini = [file("${params.hdf5_dir}/${params.hdf5_prefix}assays.h5"), file("${params.hdf5_dir}/${params.hdf5_prefix}se.rds")]
ch_hdf5 = Channel.of("${params.hdf5_prefix}").concat(Channel.fromPath(hdf5_ini))
probes = file("${params.input_probes}")
model = file("${params.trained_model}")
ch_network1 = file("${params.network1}")

include { R_SELECT_INPUTCPG_HDF5 } from '../modules/local/r_select_inputcpgs_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/gse84727/${date}"])
include { PY_RESHAPE_HDF5 } from '../modules/local/py_reshape_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/gse84727/${date}"])
include { PY_DIVIDE_TRAIN_TEST } from '../modules/local/py_divide_training_test/main.nf' addParams( options: [publish_files : ['s':'']])
include { R_WRITE_TCGA_LABELS } from '../modules/local/r_write_tcga_labels/main.nf' addParams( options: [publish_files : ['s':'']])
include { TRANSFER_LEARNING } from '../modules/local/transfer_learning/main.nf' addParams( options: [publish_dir: "gse84727/${date}"])
include { TRANSFER_LEARNING_DENSE } from '../modules/local/transfer_learning_dense/main.nf' addParams( options: [publish_dir: "gse84727_dense/${date}"])

include { PY_EXPORT_RESULTS } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "gse84727/${date}"])
include { PY_EXPORT_RESULTS  as PY_EXPORT_RESULTS_DENSE } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "gse84727_dense/${date}"])

include { R_FILTER_INVARIANTCPG_HDF5 } from '../modules/local/r_filter_invariantcpg_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/gse84727/${date}"])
include { R_SUBSTITUE_MISSING_HDF5 } from '../modules/local/r_substitue_missing_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/gse84727/${date}"])
include { TRAIN_TCGA } from '../modules/local/train_tcga/main.nf' addParams( options: [publish_dir: "gse84727/${date}"])
include { RUN_MODEL } from '../modules/local/run_model/main.nf' addParams( options: [publish_dir: "gse84727_model/${date}"])

workflow  {

  R_SELECT_INPUTCPG_HDF5( ch_hdf5.collect(), probes )
  R_WRITE_TCGA_LABELS( R_SELECT_INPUTCPG_HDF5.out.res )

  PY_RESHAPE_HDF5( R_SELECT_INPUTCPG_HDF5.out.h5, R_WRITE_TCGA_LABELS.out )

  RUN_MODEL( model, PY_RESHAPE_HDF5.out )

  PY_DIVIDE_TRAIN_TEST ( PY_RESHAPE_HDF5.out, 0.2 )
  TRANSFER_LEARNING( PY_RESHAPE_HDF5.out, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, 'transfer' )
  PY_EXPORT_RESULTS( PY_RESHAPE_HDF5.out, TRANSFER_LEARNING.out.history, TRANSFER_LEARNING.out.model, PY_DIVIDE_TRAIN_TEST.out.test )

  TRANSFER_LEARNING_DENSE( PY_RESHAPE_HDF5.out, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, model, 'transfer_dense' )
  PY_EXPORT_RESULTS_DENSE( PY_RESHAPE_HDF5.out, TRANSFER_LEARNING_DENSE.out.history, TRANSFER_LEARNING_DENSE.out.model, PY_DIVIDE_TRAIN_TEST.out.test )


}
