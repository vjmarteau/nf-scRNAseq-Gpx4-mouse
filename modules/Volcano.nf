nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process VOLCANO {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        path("*.png")

	script:
	"""
    Draw_volcano.py \\
    --adata=${adata}
    
	"""
}
