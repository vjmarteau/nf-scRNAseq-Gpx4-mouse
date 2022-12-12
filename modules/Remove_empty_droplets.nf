nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process EMPTY_DROPLETS {
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
        path(adata)

    output:
        path("cleaned_adata.h5ad"), emit: cleaned_adata

	script:
	"""
    Remove_empty_droplets.py \\
    --adata=${adata}
    
	"""
}
