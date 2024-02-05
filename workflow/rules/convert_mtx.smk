def kallisto_to_10x_naming(base_dir, filter_status):
	out = []
	if filter_status == 'filtered_gene_bc_matrices':
		subfolder = '/counts_filtered/'
	else:
		subfolder = '/counts_unfiltered/'
	out.append(base_dir + subfolder + 'cells_x_genes.mtx')
	out.append(base_dir + subfolder + 'cells_x_genes.genes.txt')
	out.append(base_dir + subfolder + 'cells_x_genes.barcodes.txt')
	out.append(base_dir + subfolder + 'cells_x_genes.genes.names.txt')
	return(out)
	
rule convert_mtx:
	input:
		lambda wildcards: kallisto_to_10x_naming(quant_path + '/quant/{SRS}/{reference}/{workflow}', wildcards.filter_status)
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/{filter_status}/matrix.mtx'
	conda:
		"../envs/convert_mtx.yaml"
	params:
		out_dir = quant_path + '/quant/{SRS}/{reference}/{workflow}/{filter_status}'
	shell:
		"""
		Rscript {git_dir}/workflow/scripts/convert_mtx.R {input} {params.out_dir}
		"""
