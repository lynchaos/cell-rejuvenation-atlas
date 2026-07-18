process DOWNLOAD_MODULE3 {
    tag 'browder2022-tms2020'
    publishDir "${params.datadir}/module3", mode: 'copy'
    label 'low'

    output:
    path 'module3'

    script:
    """
    python -m src.module3_multiomics_integration.download_data --outdir module3
    """
}

process SCVI_INTEGRATE {
    tag 'scvi-integration'
    publishDir "${params.outdir}/module3", mode: 'copy'
    label 'gpu'

    input:
    path data_dir

    output:
    path 'module3_out'

    script:
    """
    python -m src.module3_multiomics_integration.scvi_integration \
        --tms ${data_dir}/tabula_muris_senis/tms.h5ad \
        --browder ${data_dir}/browder_long7m/browder.h5ad \
        --outdir module3_out --max-epochs ${params.max_epochs}
    """
}
