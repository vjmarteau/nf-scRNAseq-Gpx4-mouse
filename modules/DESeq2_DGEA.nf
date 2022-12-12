nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process DESeq2_DGEA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        tuple val(id), path(metadata), path(count_mat)

    output:
        path("${id}.tsv"), emit: DGEA

	script:
	"""
    DESeq2_DGEA.R \\
    --count_mat=${count_mat} \\
    --metadata=${metadata} \\
    --output_file="${id}.tsv"
	"""
}
