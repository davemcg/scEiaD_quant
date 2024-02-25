args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(decontX)

input_file <- args[1]
output_file <- args[2]
sc_obj <- readRDS(input_file)

#sc_filtered_obj <- sc_obj[colData(sc_obj)$dropkick_label == 'True']
sc_filtered_obj <- subset(sc_obj, dropkick_label == 'True')

decontaminated <- decontX(as.SingleCellExperiment(sc_filtered_obj),
							background = as.SingleCellExperiment(sc_obj))

decontaminatedS <- as.Seurat(decontaminated)
saveRDS(decontaminatedS, file = output_file)
