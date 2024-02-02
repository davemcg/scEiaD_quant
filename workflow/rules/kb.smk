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
		idx = 'references/kallisto_idx/{reference}/{workflow}/index.idx',
		t2g = 'references/t2g/{reference}/{workflow}/t2g.txt',
		cdna = 'references/t2g/{reference}/{workflow}/cdna.fasta'
	conda:
		"../envs/kb.yaml"
	shell:
		"""
		kb ref --workflow {wildcards.workflow} \
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

rule kb_count:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{reference}/{workflow}/index.idx'
	output:
		h5ad = quant_path + '/quant/{SRS}/{reference}/{workflow}/counts_unfiltered/adata.h5ad',
	threads: 8
	conda:
		"../envs/kb.yaml"
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		paired_flag = lambda wildcards: '' if SRS_dict[wildcards.SRS]['umi'] else SRS_dict[wildcards.SRS]['paired'],
		out_dir = lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.reference}/{wildcards.workflow}'
	shell:
		'''
		kb count {params.paired_flag} \
					--tmp tmp{params.out_dir} \
					--workflow {wildcards.workflow} \
					-t {threads} \
					-x {params.tech} \
					-i {input.idx} \
					-o {params.out_dir} \
					--h5ad --cellranger --filter \
					{input.fastq}
		'''
