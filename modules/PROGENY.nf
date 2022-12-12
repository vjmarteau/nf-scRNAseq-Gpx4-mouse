nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process PROGENY {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)
        path(progeny)
        path(dorothea)

    output:
        path("*.png")

	script:
	"""
    PROGENy.py \\
    --adata=${adata} \\
    --progeny=${progeny} \\
    --dorothea=${dorothea}
	"""
}
