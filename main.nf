#!/usr/bin/env nextflow
/*
 * Cell Rejuvenation Atlas — main pipeline
 * Usage:
 *   nextflow run main.nf -profile test,docker                 # CI smoke test
 *   nextflow run main.nf -profile docker --module all          # full local run
 *   nextflow run main.nf -profile awsbatch --module all        # AWS Batch
 */
nextflow.enable.dsl = 2

params.module       = 'all'          // all | rejuvenation_clock | trajectory | integration | spatial | sasp
params.outdir       = 'results'
params.datadir      = 'data'
params.from_fastq   = false
params.clock_coef   = null           // path to published clock coefficients CSV
params.max_epochs   = 400

include { DOWNLOAD_MODULE1; ANALYZE_CLOCK }        from './modules/rejuvenation_clock'
include { DOWNLOAD_MODULE2; PREPROCESS; FATE }     from './modules/trajectory'
include { DOWNLOAD_MODULE3; SCVI_INTEGRATE }       from './modules/integration'
include { DOWNLOAD_MODULE4; SPATIAL }              from './modules/spatial'
include { DOWNLOAD_MODULE5; SASP }                 from './modules/sasp'

def run_module(String name) { params.module == 'all' || params.module == name }

workflow {
    if (run_module('rejuvenation_clock')) {
        m1 = DOWNLOAD_MODULE1()
        coef_ch = Channel.fromPath(params.clock_coef ?: "${projectDir}/assets/NO_COEF")
        ANALYZE_CLOCK(m1, coef_ch)
    }
    if (run_module('trajectory')) {
        m2 = DOWNLOAD_MODULE2()
        adata = PREPROCESS(m2)
        FATE(adata)
    }
    if (run_module('integration')) {
        m3 = DOWNLOAD_MODULE3()
        SCVI_INTEGRATE(m3)
    }
    if (run_module('spatial')) {
        m4 = DOWNLOAD_MODULE4()
        SPATIAL(m4)
    }
    if (run_module('sasp')) {
        m5 = DOWNLOAD_MODULE5()
        SASP(m5)
    }
}

