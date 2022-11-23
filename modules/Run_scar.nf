nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process RUN_SCAR {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(raw_adata)
        path(adata)

    output:
        path("denoised_adata.h5ad"), emit: sample
        path("filtered_adata.h5ad"), emit: filtered_adata

	script:
	"""
    Run_scar.py \\
    --raw_adata=${raw_adata} \\
    --adata=${adata}
    
	"""
}
