# Purpose

Snakemake codebase for quantifying scRNA-seq datasets for downstream plae.nei.nih.gov / scEiaD use.

Code pulled from github.com/davemcg/scEiaD 

Goal is to split apart the [scEiaD](github.com/davemcg/scEiaD) codebase into smaller, simpler units to faciliate future sample additions. The existing codebase is a bit...complicated.

# Rough Pipeline

kb count (quant)

	cellbender (empty droplet ID and ambient RNA removal)
	seurat / scanpy basic QC filtering (mito, min features, min genes)
	solo (doublet removal)
