// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_TRAIN_DNN_BIOLOGICALMODEL {

    label 'memory_medium'
    label 'cpus_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3'

    input:
    path('assay_reshaped.h5')
    path('train.pb')
    path('test.pb')
    path('params.py')
    path('inputcpgs.txt')
    path('cpgs_map.txt')
    val(name)


    output:
    path("*history_model.pb"), emit: history
    path("*labels.pb"), emit: labels
    path(name), emit: model
    path("test_list.pb"), emit: test
    path("input_cpgs.pb"), emit: input_cpgs


    script:
    """
    train_dnn_biologicalmodel.py $name $task.cpus
    """
}
