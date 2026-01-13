# Purpose

Snakemake codebase for quantifying scRNA-seq datasets for downstream [plae.nei.nih.gov](plae.nei.nih.gov) / scEiaD / [scGeyser](https://github.com/davemcg/scGeyser) use.

Broke apart the [scEiaD](github.com/davemcg/scEiaD) codebase into smaller, simpler units to faciliate future sample additions.

This repo (davemcg/scEiaD_quant) handles the quantitation (fastq -> gene x cell count matrix).

[scEiaD_modeling](https://github.com/davemcg/scEiaD_modeling) is where the species level scVI / scANVI models are built (and existing cell type labels are transferred to label-less cells) and then I (manually) iterate through the output to assign the cell type labels.

# Rough Pipeline

  - kb count (quant)
  - cellbender (empty droplet ID and ambient RNA removal)
  - seurat / scanpy basic QC filtering (mito, min features, min genes)
  - solo (doublet removal)

# Working Directories

## Run Directory
```
 /data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_2024_01

[mcgaugheyd@biowulf scEiaD_2024_01]$ cat snakejob_2025_03_20.sh
#!/bin/bash
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
mamba activate base

bash /home/mcgaugheyd/git/scEiaD_quant/Snakemake.wrapper.sh /home/mcgaugheyd/git/scEiaD_quant/workflow/Snakefile /home/mcgaugheyd/git/scEiaD_quant/config/config_2025_03_20.yaml /home/mcgaugheyd/git/scEiaD_quant/config/cluster.json
```

## Data output example (for sample SRS7483325)

```
ll /data/OGVFB_BG/scEiaD/2024_02_28/quant/SRS7483325/hs111/standard_cell/final
total 52M
drwxr-x---. 2 mcgaugheyd mcgaugheyd 4.0K Mar 20  2025 _BPdir
-rw-r-----. 1 mcgaugheyd mcgaugheyd 2.1M Mar 20  2025 noQC.seuratV5.Rdata
-rw-r-----. 1 mcgaugheyd mcgaugheyd 1.2M Mar 20  2025 QC.seuratV5.Rdata
-rw-r-----. 1 mcgaugheyd mcgaugheyd  41M Mar 20  2025 noQC.adata.h5ad
-rw-r-----. 1 mcgaugheyd mcgaugheyd 8.4M Mar 20  2025 QC.adata.solo.h5ad
-rw-r-----. 1 mcgaugheyd mcgaugheyd  590 Mar 20  2025 QC.adata.solo.obs.csv.gz
```


# Process for adding new samples

A yaml file contains the required information to run the Snakemake pipeline
```
[mcgaugheyd@biowulf scEiaD_2024_01]$ cat /home/mcgaugheyd/git/scEiaD_quant/config/config_2025_03_20.yaml
srr_sample_file: '/home/mcgaugheyd/git/scEiaD_quant/sample_meta.scEiaD_v1.2025_03_20.01.tsv.gz'
ref_file: '/home/mcgaugheyd/git/scEiaD_quant/reference_table.tsv'
git_dir: '/home/mcgaugheyd/git/scEiaD_quant/'
working_dir: '/data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_2024_01/'
fastq_path: '/data/OGVFB_BG/scEiaD'
quant_path: '/data/OGVFB_BG/scEiaD/2024_02_28'
make_10x_style_mtx: True
decontamination: "cellbender"
```

You will have to update the `srr_sample_file` with the new sample metadata
```
| sample_accession | run_accession | library_layout | reference | kb_tech |  umi | workflow | kb_sum | organism     | platform | study_accession | tissue              | sub_tissue | covariate | perturbation | inte>
| ---------------- | ------------- | -------------- | --------- | ------- | ---- | -------- | ------ | ------------ | -------- | --------------- | ------------------- | ---------- | --------- | ------------ | ---->
| SRX14524742      | SRR18390614   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Ciliary body        |            | Pt2       |              |     >
| SRX14524741      | SRR18390615   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Lens                |            | Pt14      |              |     >
| SRX14524740      | SRR18390616   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Cornea              |            | Pt14      |              |     >
| SRX14524739      | SRR18390617   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Trabecular meshwork |            | Hu0235    |              |     >
| SRX14524738      | SRR18390618   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Cornea sclera wedge |            | Hu0235    |              |     >
| SRX14524737      | SRR18390619   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Iris                |            | Hu0235    |              |     >
| SRX14524736      | SRR18390620   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Iris                |            | Hu0235    |              |     >
| SRX14524735      | SRR18390621   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Cornea              |            | Hu0235    |              |     >
| SRX14524734      | SRR18390622   | PAIRED         | hs111     | 10xv3   | True | nac      | total  | Homo sapiens | 10xv3    | SRP364915       | Ciliary body        |            | Hu0235    |              |
```

All columns (as of 2026 01 13)
```
zcat /home/mcgaugheyd/git/scEiaD_quant/sample_meta.scEiaD_v1.2025_03_20.01.tsv.gz | head -n 1
sample_accession	run_accession	library_layout	reference	kb_tech	umi	workflow	kb_sum	organism	platform	study_accession	tissue	sub_tissue	covariate	perturbation	integration_group	tissuenote	source	bam10x	comment	biosample	organ	sex	biosample_title	strain	batch	age	capture_type	enriched_cell_type	suspension_enrichment_factors	ethnicity
```

Most of these fields can be found in the [SRA run selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP364915&o=acc_s%3Aa).

Some notes for sample adding:

  1. If there are multiple run accession for a sample, then *each* run accession gets its own line
  2. the `kb_tech` column can also take a custom field (`10xv2 -w references/3M-february-2018.txt`) to give a custom white list.
  3. You can specify whether the quantitation strategy in the `kb_sum` (e.g. whether the data is single cell ("cell") or nucleus ("total").

# Running
```
cd /data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_2024_01
# see contents below.
# obviously update with the new config yaml file
sbatch snakejob_2025_03_20.sh
```

```
cat  snakejob_2025_03_20.sh
#!/bin/bash
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
mamba activate base

bash /home/mcgaugheyd/git/scEiaD_quant/Snakemake.wrapper.sh /home/mcgaugheyd/git/scEiaD_quant/workflow/Snakefile /home/mcgaugheyd/git/scEiaD_quant/config/config_2025_03_20.yaml /home/mcgaugheyd/git/scEiaD_quant/config/cluster.json
```
