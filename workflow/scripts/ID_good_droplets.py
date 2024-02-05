import scanpy as sc
import dropkick as dk

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-h", "--h5ad", help="h5ad filepath")
parser.add_argument("-o", "--csvout", help="csv output filepath")
parser.add_argument("-i", "--pngout", help="png output filepath")
args = parser.parse_args()


adata = sc.read_h5ad(args.h5ad)

adata_model = dk.dropkick(adata, n_jobs=8)
adata_filter = adata[adata.obs.dropkick_label=="True"]

qc_plt = dk.qc_summary(adata)
qc_plt.savefig(args.pngout)

obs = adata_filter.obs
obs.dropkick_label.to_csv(args.csvout)
