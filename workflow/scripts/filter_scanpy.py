import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(snakemake.input[0])
gene_tab = pd.read_csv(snakemake.input[1])

gene_tab.index = adata.var.index
adata.var['id_name'] = gene_tab['nname']

adata.var['mt'] = adata.var['id_name'].str.startswith(('MT-', 'mt-'))# annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata.write_h5ad(snakemake.output[0])

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_cells(adata, min_counts=500)


adata = adata[adata.obs.pct_counts_mt < 10, :]

adata.write_h5ad(snakemake.output[1])
