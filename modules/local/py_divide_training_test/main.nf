// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_DIVIDE_TRAIN_TEST {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3'

    input:
    path('assays.h5')
    val(prop)

    output:
    path("train.pb"), emit: train
    path("test.pb"), emit: test

    script:
    """
    divide_train_test.py $prop
    """
}
