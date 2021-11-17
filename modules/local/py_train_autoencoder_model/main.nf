// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRAIN_AUTOENCODER_MODEL {

    label 'medium_memory'
    label 'high_cpus'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.4'

    input:
    path('train.pb')
    path('test.pb')
    path('model/')
    val(name)

    output:
    path("*autoencode_history_model.pb"), emit: old_history
    path("*autoencoder_history_model.pb"), emit: auto_history
    path("*_with_autoencode"), emit: model_old
    path("*_autoencoder"), emit: model_auto

    script:
    """
    train_autoencoder.py $name $task.cpus
    """
}
