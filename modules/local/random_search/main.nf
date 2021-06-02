// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RANDOM_SEARCH {

    tag "$round"

    label 'process_high'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.0'

    input:
    path('input.pb')
    tuple val(round), path('randomconfig.py')


    output:
    path("*.pb"), emit: pickle
    path("*.tsv"), emit: tsv

    script:
    """
    randomSearchCV.py randomconfig.py $round
    """
}
