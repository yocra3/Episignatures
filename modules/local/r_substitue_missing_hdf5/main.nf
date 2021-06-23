// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process R_SUBSTITUE_MISSING_HDF5 {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignatures_rsession:1.0'

    input:
    tuple val(prefix), path(hdf5), path(rds)

    output:
    tuple val("${prefix}missingSub_"), path("*.h5"), path("*.rds"), emit: res
    path("*.h5"), emit: h5

    script:
    """
    substitute_missings.R $prefix
    """
}