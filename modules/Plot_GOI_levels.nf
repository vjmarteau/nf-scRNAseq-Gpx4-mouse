nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Plot_GOI_Levels {
    publishDir "${out_dir}", mode: "$mode"

    input:
        tuple val(id), path(metadata), path(count_mat)
        path(GOI)

    output:
        path("${id}.pdf")

	script:
	"""
    Plot_GOI_levels.R \\
    --count_mat=${count_mat} \\
    --metadata=${metadata} \\
    --GOI=${GOI} \\
    --output_file="${id}.pdf"
	"""
}
