# remove empty droplets and ambient correct counts
rule cellbender:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/adata.h5ad'
	output:
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5'),
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender.h5'),
	params:
		h5 = quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender.h5',
		working_dir = quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/',
	conda:
		"../envs/cellbender.yaml"
	shell:
		"""
		cd {params.working_dir}	
		cellbender remove-background \
			--input {input} --output {params.h5} --cuda 
		rm -f {params.working_dir}/*posterior.h5
		"""

rule h5_to_h5ad:
	input:
		quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5'
	output:
		temp(quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.h5ad'),
		quant_path + '/quant/{SRS}/{reference}/{workflow}/cellbender/cellbender_filtered.obs.csv.gz'
	conda:
		"../envs/cellbender.yaml"
	script:
		"../scripts/cellbender_h5_to_h5ad.py"

