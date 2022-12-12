nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process RUN_SCAR {
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
        path(raw_adata)
        path(adata)
        path(cell_cycle_genes)

    output:
        path("denoised_adata.h5ad")
        path("filtered_adata.h5ad"), emit: filtered_adata

	script:
	"""
    Run_scar.py \\
    --raw_adata=${raw_adata} \\
    --adata=${adata} \\
    --cell_cycle_genes=${cell_cycle_genes}
    
	"""
}
