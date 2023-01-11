nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process DENOISED_PSEUDOBULK {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        path("*samplesheet.csv"), emit: denoised_samplesheet
        path("*counts.csv"), emit: denoised_counts

	script:
	"""
    Denoised_pseudobulk.py \\
    --adata=${adata}
	"""
}
