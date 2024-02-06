import scvi 
import scanpy as sc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train(early_stopping = True, accelerator = 'gpu')

solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()

preds = solo.predict()

preds['prediction'] = solo.predict(soft = False)

preds['solo_dif'] = preds.doublet - preds.singlet

doublets = preds[(preds.prediction == 'doublet') & (preds.solo_dif > 1)]

adata.obs['solo_doublet'] = adata.obs.index.isin(doublets.index)

adata.obs['solo_score'] = preds['solo_dif']

adata.write(args.output)
