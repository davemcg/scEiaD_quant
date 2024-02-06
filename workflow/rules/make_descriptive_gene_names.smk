rule make_descriptive_gene_names:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cells_x_genes.genes.txt',
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cells_x_genes.genes.names.txt'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/gene_tab.csv'
	conda:
		'../envs/snakequant_seurat_qc.yaml'
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/make_descriptive_gene_names.R {input} {output}
		"""
