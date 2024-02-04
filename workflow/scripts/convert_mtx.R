# kallisto generates a gene x cell mtx
# packages which expect a 10X style mtx cannot handle this
# 10x style mtx is cell x gene
# so we have to import and transpose
# adapted lightly from https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/reformat.R
library(Matrix, quietly=T) 
library(DropletUtils, quietly=T)

# mtx
raw_mtx <- as(t(readMM(snakemake@input[[1]], 'CsparseMatrix') # load mtx and transpose it
# genes
rownames(raw_mtx) <- read.csv(snakemake@input[[2]], sep = '\t', header = F)[,1] # attach genes
# barcodes
colnames(raw_mtx) <- read.csv(snakemake@input[[3]], header = F, sep = '\t')[,1] # attach barcodes
# kb count 0.28 returns the names as a separate file
gene_sym <- read.csv(snakemake@input[[4]], sep = '\t', header = F)[,1]

# output dir
write10xCounts(snakemake@input[[5]], gene.symbol = gene_sym, raw_mtx, overwrite = T) # write results

