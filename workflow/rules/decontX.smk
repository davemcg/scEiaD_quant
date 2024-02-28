rule dropkick:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_unfiltered/adata.h5ad'
	output:
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick.h5ad')
	params:
		out_dir = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/'
	conda:
		"../envs/dropkick.yaml"
	shell:
		"""
		dropkick run {input} -o {params.out_dir}
		"""

rule h5ad_to_seurat:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick.h5ad'
	output:
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick.seurat.rds')
	conda:
		"../envs/sceasy.yaml"
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/convert_format.R {input} anndata {output} seurat
		"""
		
rule decontX:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick.seurat.rds'
	output:
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.seurat.rds')
	conda:
		"../envs/snakequant_seurat_qc.yaml"
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/decontX.R {input} {output}
		"""

rule seurat_to_h5ad:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.seurat.rds'	
	output:
		h5ad = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.h5ad',
		obs = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/counts_filtered/adata_dropkick_decontX.obs.csv.gz'
	conda:
		"../envs/sceasy.yaml"
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/convert_format.R {input} seurat {output.h5ad} anndata {output.obs}
		"""
