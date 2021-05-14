// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RANDOM_SEARCH {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_python:1.0'

    input:
    path('input.pb')
    path('randomconfig.py')
    val(round)

    output:
    path("*.pb"), emit: pickle
    path("*.tsv"), emit: tsv

    script:
    """
    randomSearchCV.py randomconfig.py $round
    """
}
