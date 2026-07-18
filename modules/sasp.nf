process DOWNLOAD_MODULE5 {
    tag 'basisty2020-sasp'
    publishDir "${params.datadir}/module5", mode: 'copy'
    label 'low'

    output:
    path 'module5'

    script:
    """
    python -m src.module5_proteomics_sasp.download_data --outdir module5
    """
}

process SASP {
    tag 'sasp-reanalysis'
    publishDir "${params.outdir}/module5", mode: 'copy'
    label 'medium'

    input:
    path data_dir

    output:
    path 'module5_out'

    script:
    """
    python -m src.module5_proteomics_sasp.sasp_analysis \
        --input ${data_dir} --outdir module5_out \
        ${params.sasp_design ? "--design ${params.sasp_design}" : ''}
    """
}
