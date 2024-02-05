# remove empty droplets and ambient correct counts
rule cellbender:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/adata.h5ad'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cellbender_filtered.h5'
	params:
		h5 = quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cellbender.h5',
		working_dir = quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/',
	conda:
		"../envs/cellbender.yaml"
	shell:
		"""
		cd {params.working_dir}	
		cellbender remove-background \
			--input {input} --output {params.h5} --cuda 
		"""

rule h5_to_h5ad:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cellbender_filtered.h5'
	output:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/cellbender_filtered.h5ad'
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/cellbender_h5_to_h5ad.py"
