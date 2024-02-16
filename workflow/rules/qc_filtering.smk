rule filter_seurat:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5ad',
		quant_path +'/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.obs.csv.gz',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/gene_tab.csv',
	output:
		noQC = quant_path + '/quant/{SRS}/{reference}/{workflow}/final/noQC.seuratV5.Rdata',
		QC = quant_path + '/quant/{SRS}/{reference}/{workflow}/final/QC.seuratV5.Rdata',
		faux = temp(quant_path + '/quant/{SRS}/{reference}/{workflow}/final/SEURAT_DONE')
	params:
		path = quant_path + '/quant/{SRS}/{reference}/{workflow}/final/'
	conda:
		"../envs/snakequant_seurat_qc.yaml"
	shell:
		"""
		cd {params.path}
		Rscript {git_dir}/workflow/scripts/filter_seurat.R {input} {wildcards.SRS} {output}
		cd -
		touch {output.faux}
		"""

rule filter_scanpy:
	input:
		quant_path +'/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/gene_tab.csv',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/SEURAT_DONE'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/noQC.adata.h5ad',
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}/final/QC.adata.h5ad')
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/filter_scanpy.py"
		
