#!/bin/bash
# Run scRNAseq gpx4 analysis pipeline

nextflow run ./main.nf -profile cluster \
-resume