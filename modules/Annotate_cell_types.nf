nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process ANNOTATE_CELL_TYPES {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)
        path(marker_genes)

    output:
        path("annotated_adata.h5ad"), emit: annotated_adata
        path("*samplesheet.csv"), emit: samplesheet
        path("*counts.csv"), emit: counts

	script:
	"""
    Annotate_cell_types.py \\
    --adata=${adata} \\
    --marker_genes=${marker_genes}
	"""
}
