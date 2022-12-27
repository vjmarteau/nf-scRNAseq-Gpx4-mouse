nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process CONVERT_SEURAT {
    publishDir "${out_dir}", mode: "$mode"

    input:
         path(counts_matrix)
         path(denoised_matrix)
         path(features)
         path(barcodes)
         path(metadata)
         path(metadata_var)
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
