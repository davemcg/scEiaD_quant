library(dplyr)
library(BPCells)
library(Seurat) #version >= 5
args = commandArgs(trailingOnly=TRUE)

adata <- open_matrix_anndata_hdf5(args[1])

bp_out <- gsub(args[1], ".h5ad", "_BPdir")
write_matrix_dir(mat = adata, dir = bp_out, overwrite =TRUE)
   
mat <- open_matrix_dir(dir = bp_out)

metadata <- read.csv(args[2])

gene_tab <- read.csv(args[3], sep = ',')

row.names(mat) <- gene_tab$nname
seurat <- CreateSeuratObject(counts = list(mat), meta.data = metadata)

seurat@meta.data$orig.ident = args[4]

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-|^mt-")

mincount <- 500
minfeature <- 300
mito <- 10
seurat_filter <- subset(seurat, subset=nCount_RNA>=mincount & nFeature_RNA>=minfeature & percent.mt<=mito)

save(seurat, file = args[5])
save(seurat_filter, file = args[6])

