args = commandArgs(trailingOnly=TRUE)

library(sceasy)
input_file = args[1]
input_format = args[2]
output_file = args[3]
output_format = args[4]

print(args)
if (input_format == 'anndata'){
	sceasy::convertFormat(input_file, from="anndata", to=output_format,
                       outFile=output_file)
} else {
	input_obj <- readRDS(input_file)
	sceasy::convertFormat(input_obj, from=input_format, to=output_format,
                       outFile=output_file)
	if (!is.null(args[5])){
		print("Write OBS")
		if (input_format == 'seurat'){
			meta = input_obj@meta.data
			write.csv(meta, file = gzfile(args[5]))
		} else if (input_format == 'sce'){
			meta = colData(input_obj)
			write.csv(meta, file = gzfile(args[5]))
		}
	}
}
