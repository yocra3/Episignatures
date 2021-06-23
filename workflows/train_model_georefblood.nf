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
ch_random2 = file("${params.config2}")
ch_random3 = file("${params.config3}")

ch_network1 = file("${params.network1}")

ch_round1 = Channel.of('round1').concat(Channel.fromPath(ch_random1))
ch_round2 = Channel.of('round2').concat(Channel.fromPath(ch_random2))
ch_round3 = Channel.of('round3').concat(Channel.fromPath(ch_random3))

//ch_random = ch_round1.collect().concat(ch_round2.collect()).concat(ch_round3.collect())
//ch_random = ch_round1.collect().concat(ch_round2.collect().ifEmpty([]))

include { PY_DIVIDE_TRAIN_TEST } from '../modules/local/py_divide_training_test/main.nf' addParams( options: [publish_files : ['s':'']])
include { RANDOM_SEARCH_2DNN } from '../modules/local/random_search_twodnn/main.nf' addParams( options: [publish_dir: "GEOrefblood_model/${date}"])
include { TRAIN_TCGA } from '../modules/local/train_tcga/main.nf' addParams( options: [publish_dir: "GEOrefblood_model/${date}"])
include { PY_EXPORT_RESULTS } from '../modules/local/py_export_results/main.nf' addParams( options: [publish_dir: "GEOrefblood_model/${date}"])

workflow  {

  PY_DIVIDE_TRAIN_TEST ( ch_hdf5, 0.2 )

  RANDOM_SEARCH_2DNN( PY_DIVIDE_TRAIN_TEST.out.train.collect(), ch_round1.collect())

//  TRAIN_TCGA( ch_hdf5, PY_DIVIDE_TRAIN_TEST.out.train, PY_DIVIDE_TRAIN_TEST.out.test, ch_network1, 'model3' )
//  PY_EXPORT_RESULTS( ch_hdf5, TRAIN_TCGA.out.history, TRAIN_TCGA.out.model, PY_DIVIDE_TRAIN_TEST.out.test )

}
