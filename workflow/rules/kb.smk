rule download_ref:
	output:
		gtf = 'references/gtf/{reference}.gtf.gz',
		fasta = 'references/fasta/{reference}.fa.gz'
	params:
		gtf_url = lambda wildcards: ref_dict[wildcards.reference]['gtf'],
		fasta_url = lambda wildcards: ref_dict[wildcards.reference]['fasta']
	shell:
		'''
		mkdir -p references/gtf
		mkdir -p references/fasta
		wget -O {output.gtf} {params.gtf_url}
		wget -O {output.fasta} {params.fasta_url}
		'''

rule kb_ref:
	input:
		'references/fasta/{reference}.fa.gz',
		'references/gtf/{reference}.gtf.gz'
	output:
		idx = 'references/kallisto_idx/{reference}/standard/index.idx',
		t2g = 'references/t2g/{reference}/standard/t2g.txt',
		cdna = 'references/t2g/{reference}/standard/cdna.fasta'
	conda:
		"../envs/kb.yaml"
	shell:
		"""
		rm -fr tmp{wildcards.reference}{wildcards.workflow}
		kb ref --workflow {wildcards.workflow} \
			--tmp tmp{wildcards.reference}{wildcards.workflow} \
			-i {output.idx} \
			-g {output.t2g} \
			-f1 {output.cdna} \
				--include-attribute gene_biotype:protein_coding \
				--include-attribute gene_biotype:lncRNA \
  				--include-attribute gene_biotype:lincRNA \
  				--include-attribute gene_biotype:antisense \
  				--include-attribute gene_biotype:IG_LV_gene \
  				--include-attribute gene_biotype:IG_V_gene \
  				--include-attribute gene_biotype:IG_V_pseudogene \
  				--include-attribute gene_biotype:IG_D_gene \
  				--include-attribute gene_biotype:IG_J_gene \
  				--include-attribute gene_biotype:IG_J_pseudogene \
  				--include-attribute gene_biotype:IG_C_gene \
  				--include-attribute gene_biotype:IG_C_pseudogene \
  				--include-attribute gene_biotype:TR_V_gene \
  				--include-attribute gene_biotype:TR_V_pseudogene \
  				--include-attribute gene_biotype:TR_D_gene \
  				--include-attribute gene_biotype:TR_J_gene \
  				--include-attribute gene_biotype:TR_J_pseudogene \
  				--include-attribute gene_biotype:TR_C_gene \
			{input}
		"""

rule kb_ref_nac:
	input:
		'references/fasta/{reference}.fa.gz',
		'references/gtf/{reference}.gtf.gz'
	output:
		idx = 'references/kallisto_idx/{reference}/nac/index.idx',
		t2g = 'references/t2g/{reference}/nac/t2g.txt',
		cdna = 'references/t2g/{reference}/nac/cdna.fasta',
		f2 = 'references/t2g/{reference}/nac/unprocessed.fasta',
		c1 = 'references/t2g/{reference}/nac/t2c.cdna.txt',
		c2 = 'references/t2g/{reference}/nac/t2c.unprocessed.txt'
	conda:
		"../envs/kb.yaml"
	shell:
		"""
		rm -fr tmp{wildcards.reference}nac
		kb ref --workflow nac \
			--tmp tmp{wildcards.reference}nac \
			-i {output.idx} \
			-g {output.t2g} \
			-f1 {output.cdna} \
			-f2 {output.f2} \
			-c1 {output.c1} \
			-c2 {output.c2} \
				--include-attribute gene_biotype:protein_coding \
				--include-attribute gene_biotype:lncRNA \
  				--include-attribute gene_biotype:lincRNA \
  				--include-attribute gene_biotype:antisense \
  				--include-attribute gene_biotype:IG_LV_gene \
  				--include-attribute gene_biotype:IG_V_gene \
  				--include-attribute gene_biotype:IG_V_pseudogene \
  				--include-attribute gene_biotype:IG_D_gene \
  				--include-attribute gene_biotype:IG_J_gene \
  				--include-attribute gene_biotype:IG_J_pseudogene \
  				--include-attribute gene_biotype:IG_C_gene \
  				--include-attribute gene_biotype:IG_C_pseudogene \
  				--include-attribute gene_biotype:TR_V_gene \
  				--include-attribute gene_biotype:TR_V_pseudogene \
  				--include-attribute gene_biotype:TR_D_gene \
  				--include-attribute gene_biotype:TR_J_gene \
  				--include-attribute gene_biotype:TR_J_pseudogene \
  				--include-attribute gene_biotype:TR_C_gene \
			{input}
		"""

rule kb_umi_count:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{reference}/standard/index.idx',
		t2g = 'references/t2g/{reference}/standard/t2g.txt',
		cdna = 'references/t2g/{reference}/standard/cdna.fasta'
	output:
		quant_path + '/quant/{SRS}/{reference}/standard_{sum}/counts_unfiltered/cells_x_genes.mtx',
		quant_path + '/quant/{SRS}/{reference}/standard_{sum}/counts_unfiltered/cells_x_genes.genes.txt',
		quant_path + '/quant/{SRS}/{reference}/standard_{sum}/counts_unfiltered/cells_x_genes.barcodes.txt',
		quant_path + '/quant/{SRS}/{reference}/standard_{sum}/counts_unfiltered/cells_x_genes.genes.names.txt',
		quant_path + '/quant/{SRS}/{reference}/standard_{sum}/counts_unfiltered/adata.h5ad'
	threads: 8
	conda:
		"../envs/kb.yaml"
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		out_dir = lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.reference}/standard_{wildcards.sum}'
	shell:
		'''
		rm -fr tmp{wildcards.SRS}{wildcards.reference}standard
		kb count \
			--tmp tmp{wildcards.SRS}{wildcards.reference}standard \
			--workflow standard \
			--sum {wildcards.sum} \
			-g {input.t2g} \
			-t {threads} \
			-x {params.tech} \
			-i {input.idx} \
			-o {params.out_dir} \
			--h5ad \
			{input.fastq}
		rm -f {params.out_dir}*bus
		'''

rule kb_umi_count_nac:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{reference}/nac/index.idx',
		t2g = 'references/t2g/{reference}/nac/t2g.txt',
		cdna = 'references/t2g/{reference}/nac/cdna.fasta',
		f2 = 'references/t2g/{reference}/nac/unprocessed.fasta',
		c1 = 'references/t2g/{reference}/nac/t2c.cdna.txt',
		c2 = 'references/t2g/{reference}/nac/t2c.unprocessed.txt'
	output:
		a = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/cells_x_genes.mature.mtx',
		b = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/cells_x_genes.nascent.mtx',
		c = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/cells_x_genes.genes.txt',
		d = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/cells_x_genes.barcodes.txt',
		e = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/cells_x_genes.genes.names.txt',
		f = quant_path + '/quant/{SRS}/{reference}/nac_{sum}/counts_unfiltered/adata.h5ad'
	threads: 8
	conda:
		"../envs/kb.yaml"
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		out_dir = lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.reference}/nac_{wildcards.sum}'
	shell:
		'''
		rm -fr tmp{wildcards.SRS}{wildcards.reference}nac
		kb count \
			--tmp tmp{wildcards.SRS}{wildcards.reference}nac \
			--workflow nac \
			--sum {wildcards.sum} \
			-g {input.t2g} \
			-t {threads} \
			-x {params.tech} \
			-i {input.idx} \
			-c1 {input.c1} \
			-c2 {input.c2} \
			-o {params.out_dir} \
			--h5ad \
			{input.fastq}
		rm -f {params.out_dir}*bus
		mv {output.f} {output.f}O
		python {git_dir}/workflow/scripts/kb_h5ad_agg.py {output.f}O {output.f}
		'''


rule kb_bulk_count:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{reference}/{workflow}/index.idx',
		t2g = 'references/t2g/{reference}/{workflow}/t2g.txt',
		cdna = 'references/t2g/{reference}/{workflow}/cdna.fasta'
	output:
		ec = quant_path + '/quant/{SRS}/{reference}/{workflow}_{sum}/matrix.ec',
	threads: 1 
	conda:
		"../envs/kb.yaml"
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		paired_flag = lambda wildcards: SRS_dict[wildcards.SRS]['parity'],
		out_dir = lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.reference}/{wildcards.workflow}_{wildcards.sum}'
	shell:
		'''
		rm -fr tmp{wildcards.SRS}{wildcards.reference}{wildcards.workflow}
		kb count {params.paired_flag} \
					--matrix-to-files \
					--tmp tmp{wildcards.SRS}{wildcards.reference}{wildcards.workflow} \
					--workflow {wildcards.workflow} \
					-g {input.t2g} \
					-t {threads} \
					-x {params.tech} \
					-i {input.idx} \
					-o {params.out_dir} \
					{input.fastq}
		rm -f {params.out_dir}/*bus
		'''
