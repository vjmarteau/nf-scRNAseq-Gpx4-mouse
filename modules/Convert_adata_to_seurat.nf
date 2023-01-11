nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process CONVERT_EXTRACT_ADATA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        tuple path("counts_matrix.mtx"),
        path("denoised_matrix.mtx"),
        path("features.tsv"),
        path("barcodes.tsv"),
        path("metadata.tsv"),
        path("metadata_var.tsv"),
        path("umap.tsv"), emit: convert
    

	script:
	"""
    Convert.py \\
    --adata=${adata}
    
	"""
}

process CONVERT_ASSEMBLE_SEURAT {
    publishDir "${out_dir}", mode: "$mode"

    input:
         tuple path(counts_matrix),
         path(denoised_matrix),
         path(features),
         path(barcodes),
         path(metadata),
         path(metadata_var),
         path(umap)

    output:
        path("*.rds"), emit: seurat

	script:
	"""
    Convert.R \\
    --counts_mtx=${counts_matrix} \\
    --denoised_mtx=${denoised_matrix} \\
    --features=${features} \\
    --barcodes=${barcodes} \\
    --metadata=${metadata} \\
    --metadata_var=${metadata_var} \\
    --umap_tsv=${umap}
    
	"""
}

workflow CONVERT_ADATA_TO_SEURAT {
    take:
        adata

    main:
        CONVERT_EXTRACT_ADATA(adata)
        CONVERT_ASSEMBLE_SEURAT(CONVERT_EXTRACT_ADATA.out.convert)

    emit:
        seurat = CONVERT_ASSEMBLE_SEURAT.out.seurat
}
