nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process CONVERT_ADATA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        path("counts_matrix.mtx"), emit: counts_matrix
        path("denoised_matrix.mtx"), emit: denoised_matrix
        path("features.tsv"), emit: features
        path("barcodes.tsv"), emit: barcodes
        path("metadata.tsv"), emit: metadata
        path("metadata_var.tsv"), emit: metadata_var
        path("umap.tsv"), emit: umap
    

	script:
	"""
    Convert.py \\
    --adata=${adata}
    
	"""
}
