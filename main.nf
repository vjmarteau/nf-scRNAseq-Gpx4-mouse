#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LOAD_ADATA } from "./modules/Load_adata"
include { RUN_SCAR } from "./modules/Run_scar"
include { RUN_SCVI_AND_SOLO } from "./modules/Run_scVI_and_Solo"
include { ANNOTATE_CELL_TYPES } from "./modules/Annotate_cell_types"
include { DESeq2_DGEA } from "./modules/DESeq2_DGEA"

include { Render_py } from "./modules/Render_py"
include { Render_R } from "./modules/Render_R"


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
    
    ch_samplesheets_by_cell_type = ANNOTATE_CELL_TYPES.out.samplesheet.flatten().map { 
        it -> [it.baseName.replace("_samplesheet", ""), it]
    }.view()

    ch_counts_by_cell_type = ANNOTATE_CELL_TYPES.out.counts.flatten().map { 
        it -> [it.baseName.replace("_counts", ""), it]
    }.view()

    ch_deseq2_input = ch_samplesheets_by_cell_type.join(ch_counts_by_cell_type).view()

    DESeq2_DGEA(ch_deseq2_input)
    //DESeq2_DGEA(Reformat_data.out.all.mix(Reformat_data2.out.all.view()), prefix)

    //channel1.map{ it -> [it.baseName, it]}
    //join
    //mix
    //collect
    //flatten

    }
