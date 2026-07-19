process DOWNLOAD_MODULE4 {
    tag 'allen2023-merfish'
    publishDir "${params.datadir}/module4", mode: 'copy'
    label 'low'

    output:
    path 'module4'

    script:
    """
    python -m src.module4_spatial_aging.download_data --outdir module4
    """
}

process SPATIAL {
    tag 'merfish-spatial'
    publishDir "${params.outdir}/module4", mode: 'copy'
    label 'medium'

    input:
    path data_dir

    output:
    path 'module4_out'

    script:
    """
    python -m src.module4_spatial_aging.spatial_analysis \
        --input ${data_dir} --outdir module4_out
    """
}
