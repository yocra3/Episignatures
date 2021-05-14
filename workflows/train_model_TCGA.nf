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

ch_hdf5 = Channel.fromPath(file("${params.hdf5_file}"))
ch_random1 = file("${params.config1}")
ch_network1 = file("${params.network1}")


include { PY_DIVIDE_TRAIN_TEST } from '../modules/local/py_divide_training_test/main.nf' addParams( options: [publish_files : ['s':'']])
include { RANDOM_SEARCH } from '../modules/local/random_search/main.nf' addParams( options: [publish_dir: "tcga_model/${date}"])
include { TRAIN_TCGA } from '../modules/local/train_tcga/main.nf' addParams( options: [publish_dir: "tcga_model/${date}"])
include { PY_EXPORT_RESULTS } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "tcga_model/${date}"])

workflow  {

  PY_DIVIDE_TRAIN_TEST ( ch_hdf5 )
  RANDOM_SEARCH( PY_DIVIDE_TRAIN_TEST.out.train, ch_random1, 'round1' )
  TRAIN_TCGA( PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, ch_network1, 'round1' )
  PY_EXPORT_RESULTS( ch_hdf5, TRAIN_TCGA.out.history, TRAIN_TCGA.out.model, PY_DIVIDE_TRAIN_TEST.out.test )

}
