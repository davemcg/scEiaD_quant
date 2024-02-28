rule filter_seurat_cellbender:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/cellbender/cellbender_filtered.h5ad',
		quant_path +'/quant/{SRS}/{reference}/{workflow}_{sum}/cellbender/cellbender_filtered.obs.csv.gz',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_unfiltered/gene_tab.csv',
	output:
		noQC = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/noQC.seuratV5.Rdata',
		QC = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.seuratV5.Rdata',
		faux = temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/SEURAT_DONE')
	params:
		path = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/'
	conda:
		"../envs/snakequant_seurat_qc.yaml"
	shell:
		"""
		cd {params.path}
		Rscript {git_dir}/workflow/scripts/filter_seurat.R {input} {wildcards.SRS} {output}
		cd -
		touch {output.faux}
		"""

rule filter_scanpy_cellbender:
	input:
		quant_path +'/quant/{SRS}/{reference}/{workflow}_{sum}/cellbender/cellbender_filtered.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_unfiltered/gene_tab.csv',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/SEURAT_DONE'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/noQC.adata.h5ad',
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.adata.h5ad')
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/filter_scanpy.py"

		
rule filter_seurat_decontX:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.h5ad',
		quant_path +'/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.obs.csv.gz',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_unfiltered/gene_tab.csv',
	output:
		noQC = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/noQC.seuratV5.Rdata',
		QC = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.seuratV5.Rdata',
		faux = temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/SEURAT_DONE')
	params:
		path = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/'
	conda:
		"../envs/snakequant_seurat_qc.yaml"
	shell:
		"""
		cd {params.path}
		Rscript {git_dir}/workflow/scripts/filter_seurat.R {input} {wildcards.SRS} {output}
		rm -rf _BPdir
		cd -
		touch {output.faux}
		"""

rule filter_scanpy_decontX:
	input:
		quant_path +'/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_unfiltered/gene_tab.csv',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/SEURAT_DONE'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/noQC.adata.h5ad',
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.adata.h5ad')
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/filter_scanpy.py"
