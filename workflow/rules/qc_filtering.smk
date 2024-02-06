rule filter_seurat:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5ad',
		quant_path +'/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.obs.csv.gz',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/gene_tab.csv',
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/noQC.seuratV5.Rdata',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/QC.seuratV5.Rdata'
	conda:
		"../envs/snakequant_seurat_qc.yaml"
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/filter_seurat.R {input} {wildcards.SRS} {output}
		"""

rule filter_scanpy:
	input:
		quant_path +'/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/gene_tab.csv'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/noQC.adata.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/final/QC.adata.h5ad'
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/filter_scanpy.py"
		
