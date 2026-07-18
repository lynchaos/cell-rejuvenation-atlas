process DOWNLOAD_MODULE2 {
    tag 'schiebinger2019'
    publishDir "${params.datadir}/module2", mode: 'copy'
    label 'low'

    output:
    path 'module2'

    script:
    """
    python -m src.module2_reprogramming_trajectory.download_data --outdir module2
    """
}

process PREPROCESS {
    tag 'scanpy-preprocess'
    publishDir "${params.outdir}/module2", mode: 'copy'
    label 'highmem'

    input:
    path data_dir

    output:
    path 'adata.h5ad'

    script:
    """
    python -m src.module2_reprogramming_trajectory.preprocess \
        --input ${data_dir} --out adata.h5ad
    """
}

process FATE {
    tag 'wot-fate-analysis'
    publishDir "${params.outdir}/module2", mode: 'copy'
    label 'medium'

    input:
    path adata

    output:
    path 'module2_out'

    script:
    """
    python -m src.module2_reprogramming_trajectory.fate_analysis \
        --adata ${adata} --outdir module2_out
    """
}
