// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PY_EXTRACT_BIOMODEL_FEATURES {

    label 'memory_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.3'

    input:
    path('input_list.pb')
    path('model/')

    output:
    path("*.tsv"), emit: res

    script:
    """
    extract_features_bio.py
    """
}
