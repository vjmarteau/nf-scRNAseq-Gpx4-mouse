nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process UMAP {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)
        path(marker_genes)

    output:
        path("*.png")

	script:
	"""
    Plot_umap.py \\
    --adata=${adata} \\
    --marker_genes=${marker_genes}
    
	"""
}
