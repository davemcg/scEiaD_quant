import scanpy as sc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.X = adata.layers['ambiguous'] + adata.layers['mature'] + adata.layers['nascent']

adata.write_h5ad(args.output)
