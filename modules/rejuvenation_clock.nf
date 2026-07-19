process DOWNLOAD_MODULE1 {
    tag 'gill2022'
    publishDir "${params.datadir}/module1", mode: 'copy'
    label 'low'

    output:
    path 'module1'

    script:
    """
    python -m src.module1_rejuvenation_clock.download_data --outdir module1
    """
}

process ANALYZE_CLOCK {
    tag 'epigenetic-clock'
    publishDir "${params.outdir}/module1", mode: 'copy'
    label 'medium'

    input:
    path data_dir
    path coef

    output:
    path 'module1_out'

    script:
    def coef_arg = coef.name != 'NO_COEF' ? "--coef ${coef}" : ''
    """
    beta_arg="${data_dir}/module1_beta.csv"
    if [ ! -f "\$beta_arg" ]; then
        beta_arg="${data_dir}/GSE165179"
    fi

    train_arg=""
    if [ "${coef.name}" = "NO_COEF" ] && [ -f "${data_dir}/module1_train_ages.csv" ]; then
        train_arg="--train-beta ${data_dir}/module1_beta.csv --train-ages ${data_dir}/module1_train_ages.csv"
    fi

    python -m src.module1_rejuvenation_clock.analyze_gill \
        --beta \$beta_arg \
        --metadata ${data_dir}/metadata.csv \
        ${coef_arg} \
        \$train_arg \
        --outdir module1_out
    """
}
