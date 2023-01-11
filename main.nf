#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LOAD_ADATA } from "./modules/Load_adata"
include { EMPTY_DROPLETS } from "./modules/Remove_empty_droplets"
include { RUN_SCAR } from "./modules/Run_scar"
include { RUN_SCVI_AND_SOLO } from "./modules/Run_scVI_and_Solo"
include { ANNOTATE_CELL_TYPES } from "./modules/Annotate_cell_types"
include { DENOISED_PSEUDOBULK } from "./modules/Denoised_pseudobulk"
include { DESeq2_DGEA } from "./modules/DESeq2_DGEA"
include { DESeq2_DGEA as DESeq2_DGEA_DENOISED} from "./modules/DESeq2_DGEA"
include { VOLCANO } from "./modules/Volcano"
include { UMAP } from "./modules/Plot_umap"
include { PROGENY } from "./modules/PROGENY"
include { Plot_GOI_Levels } from "./modules/Plot_GOI_levels"
include { JUPYTERNOTEBOOK as JUPYTER_TEST } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as COMPOSITION } from "./modules/local/jupyternotebook/main"
include { RMARKDOWNNOTEBOOK as RMARKDOWN_TEST } from "./modules/local/rmarkdownnotebook/main"
include { CONVERT_ADATA_TO_SEURAT } from "./modules/Convert_adata_to_seurat"
//include { CONVERT_SEURAT as CONVERT_SEURAT2 } from "./modules/Convert_Seurat"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    cell_cycle_genes = file(params.cell_cycle_genes, checkIfExists: true)
    marker_genes = file(params.marker_genes, checkIfExists: true)
    progeny = file(params.progeny, checkIfExists: true)
    dorothea = file(params.dorothea, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)
    GOI = file(params.GOI, checkIfExists: true)

    // start workflow
    LOAD_ADATA(samplesheet, ch_input_files)
    EMPTY_DROPLETS(LOAD_ADATA.out.adata)
    RUN_SCAR(LOAD_ADATA.out.raw_adata, EMPTY_DROPLETS.out.cleaned_adata, cell_cycle_genes)
    RUN_SCVI_AND_SOLO(RUN_SCAR.out.filtered_adata)
    //ANNOTATE_CELL_TYPES2(RUN_SCVI_AND_SOLO.out.integrated_adata, marker_genes)
    ANNOTATE_CELL_TYPES(RUN_SCVI_AND_SOLO.out.integrated_adata, marker_genes)
    
    ch_samplesheets_by_cell_type = ANNOTATE_CELL_TYPES.out.samplesheet.flatten().map { 
        it -> [it.baseName.replace("_samplesheet", ""), it]
    }//.view()

    ch_counts_by_cell_type = ANNOTATE_CELL_TYPES.out.counts.flatten().map {
        it -> [it.baseName.replace("_counts", ""), it]
    }//.view()

    ch_deseq2_input = ch_samplesheets_by_cell_type.join(ch_counts_by_cell_type)//.view()
    DESeq2_DGEA(ch_deseq2_input)

    // Get denosied pseudobulk DGEA
    DENOISED_PSEUDOBULK(ANNOTATE_CELL_TYPES.out.annotated_adata)

    ch_samplesheets_by_cell_type_denoised = DENOISED_PSEUDOBULK.out.denoised_samplesheet.flatten().map { 
        it -> [it.baseName.replace("_denoised_samplesheet", ""), it]
    }//.view()

    ch_counts_by_cell_type_denoised = DENOISED_PSEUDOBULK.out.denoised_counts.flatten().map {
        it -> [it.baseName.replace("_denoised_counts", ""), it]
    }//.view()

    ch_deseq2_input_denoised = ch_samplesheets_by_cell_type_denoised.join(ch_counts_by_cell_type_denoised)//.view()
    DESeq2_DGEA_DENOISED(ch_deseq2_input_denoised)


    VOLCANO(ANNOTATE_CELL_TYPES.out.annotated_adata)

    UMAP(ANNOTATE_CELL_TYPES.out.annotated_adata, marker_genes)

    PROGENY(ANNOTATE_CELL_TYPES.out.annotated_adata, progeny, dorothea)
    
    COMPOSITION(
        Channel.value([
            [id: "02_cell_type_composition"],
            file("${projectDir}/analysis/03-cell_type_composition.py", checkIfExists: true)
        ]),
        Channel.value(
            ["adata_path": "annotated_adata.h5ad"]
        ),
        ANNOTATE_CELL_TYPES.out.annotated_adata
    )

    CONVERT_ADATA_TO_SEURAT(ANNOTATE_CELL_TYPES.out.annotated_adata)
    
    //CONVERT_SEURAT2(CONVERT_ADATA.out.convert)

     //Plot_GOI_Levels(ch_deseq2_input, GOI)

    //JUPYTER_TEST(
    //    Channel.value([
    //        [id: "01_test"],
    //        file("${projectDir}/analysis/test_render_notebook.ipynb", checkIfExists: true)
    //    ]),
    //    Channel.value(
    //        ["input_adata": "adata_nodoublet.h5ad"]
    //    ),
    //    RUN_SCVI_AND_SOLO.out.adata_nodoublet
    //)
    //RMARKDOWN_TEST(
    //    Channel.value([
    //        [id: "02_test_rmd"],
    //        file("${projectDir}/analysis/test_render_R.Rmd", checkIfExists: true)
    //    ]),
    //    Channel.value(
    //        ["samplesheet": "diff_genes_Paneth_counts.csv"]
    //    ),
    //    Channel.fromPath("$projectDir/results/scrnaseq/data/pseudo/diff_genes_Paneth_counts.csv")
    //)

    //DESeq2_DGEA(Reformat_data.out.all.mix(Reformat_data2.out.all.view()), prefix)

    //channel1.map{ it -> [it.baseName, it]}
    //join
    //mix
    //collect
    //flatten
 // RUN_SCVI_AND_SOLO.out.adata_nodoublet.mix(LOAD_ADATA.out.raw_adata).mix().collect()
    }
