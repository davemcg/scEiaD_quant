rule ID_good_droplets:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/{filter_status}/cells_x_genes.mtx'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/{filter_status}/dropkick_good_barcodes.csv.gz',
		quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/{filter_status}/dropkick_qc.png'
	threads: 8
	conda:
		"../envs/dropkick.yaml"
	shell:
		"""
		python ../scripts/dropkick.py {input} {output}
		"""
