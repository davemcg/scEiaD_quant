rule solo:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.adata.h5ad'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.adata.solo.h5ad',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/final/QC.adata.solo.obs.csv.gz'
	shell:
		"""
		module load scvitools/1.0.4.gpu; python-scvitools {git_dir}/workflow/scripts/solo.py {input} {output}
		"""
