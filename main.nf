#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LOAD_ADATA } from "./modules/Load_adata"
include { RUN_SCAR } from "./modules/Run_scar"
include { RUN_SCVI_AND_SOLO } from "./modules/Run_scVI_and_Solo"
include { ANNOTATE_CELL_TYPES } from "./modules/Annotate_cell_types"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    marker_genes = file(params.marker_genes, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)

    // start workflow
    LOAD_ADATA(samplesheet, ch_input_files)
    RUN_SCAR(LOAD_ADATA.out.raw_adata, LOAD_ADATA.out.adata)
    RUN_SCVI_AND_SOLO(RUN_SCAR.out.filtered_adata)
    ANNOTATE_CELL_TYPES(RUN_SCVI_AND_SOLO.out.adata_nodoublet2, marker_genes)
    }
