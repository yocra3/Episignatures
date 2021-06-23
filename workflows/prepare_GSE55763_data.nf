#!/usr/bin/env nextflow
/*
========================================================================================
Prepare data from GSE55763 for training
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////
hdf5_ini = [file("${params.hdf5_dir}/${params.hdf5_prefix}_assays.h5"), file("${params.hdf5_dir}/${params.hdf5_prefix}_se.rds")]
ch_hdf5 = Channel.of("${params.hdf5_prefix}_").concat(Channel.fromPath(hdf5_ini))

date = java.time.LocalDate.now()

include { R_FILTER_MISSINGS_HDF5 } from '../modules/local/r_filter_missings_hdf5/main.nf' addParams( options: [publish_files : ['s':'']])
include { R_SELECT_AUTOSOMICCPG_HDF5 } from '../modules/local/r_select_autosomic_cpg_hdf5/main.nf' addParams( options: [publish_files : ['s':'']])
include { R_SUBSTITUE_MISSING_HDF5 } from '../modules/local/r_substitue_missing_hdf5/main.nf' addParams( options: [publish_files : ['s':'']])
include { R_CPGS_MEDIANS } from '../modules/local/r_cpgs_medians/main.nf' addParams( options: [publish_dir: "preprocess/GSE55763/${date}"])
include { R_WRITE_TCGA_LABELS } from '../modules/local/r_write_tcga_labels/main.nf' addParams( options: [publish_files : ['s':'']])

include { PY_RESHAPE_HDF5 } from '../modules/local/py_reshape_hdf5/main.nf' addParams( options: [publish_dir: "preprocess/GSE55763/${date}"])


workflow  {

  R_SELECT_AUTOSOMICCPG_HDF5( ch_hdf5.collect() )
  R_FILTER_MISSINGS_HDF5( R_SELECT_AUTOSOMICCPG_HDF5.out )
  R_SUBSTITUE_MISSING_HDF5( R_FILTER_MISSINGS_HDF5.out )

  R_CPGS_MEDIANS( R_FILTER_MISSINGS_HDF5.out.res )
  R_WRITE_TCGA_LABELS( ch_hdf5.collect() )

  PY_RESHAPE_HDF5( R_SUBSTITUE_MISSING_HDF5.out.h5, R_WRITE_TCGA_LABELS.out )
}
