from cellbender.remove_background.downstream import anndata_from_h5
import scanpy as sc
from scipy.sparse import csr_matrix

adata = anndata_from_h5(snakemake.input[0])
adata.X = csr_matrix(adata.X)
adata.write_h5ad(snakemake.output[0])
