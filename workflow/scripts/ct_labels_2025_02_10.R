# fix srp510712 nrpc getting labelled as RPC 
# (should be neurogenic)

ct <- data.table::fread('/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2025_02_04.02.csv.gz')

#  ct %>% mutate(MajorCellType = case_when(SubCellType == 'NRPC' ~ 'neurogenic',
#                                       TRUE ~ MajorCellType)) %>% 
#  write_csv('/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2025_02_04.03.csv.gz')



# 2025 02 24

ct <- data.table::fread('/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2025_02_04.03.csv.gz')
smeta <- data.table::fread('sample_meta.scEiaD_v1.2025_02_03.02.tsv.gz')
# found a way to bring in the "missing" labels from srp443999
# they labelled their seurat_clusters in a separate excel sheet (sup from paper)

# srp443999
srp443999_files <- list.files('data/srp443999', full.names = TRUE)
srp443999_list <- set_names(srp443999_files) %>% 
  map(., read_csv)
srp443999_cells <- list()
for (i in names(srp443999_list)){
  srp443999_cells[[i]] <- try({srp443999_list[[i]] %>% select(`...1`, seurat_clusters, CellType)})
}

# grab biosample titles
smeta_srp443999 <- smeta %>% filter(study_accession == "SRP443999")
samn <- smeta_srp443999 %>% pull(biosample)
attrs <- attribute_df_maker(samn, delay = 1.4)
smeta_srp443999$title <- attrs %>% filter(attribute == 'biosample_title' ) %>% pull(value)

srp443999_cells_cm <- srp443999_cells %>% 
  bind_rows(.id = 'file') %>% 
  mutate(
    MajorCellType = case_when(CellType == 'ACs' ~ 'amacrine',
                              #CellType == 'Amacrine and horizontal cells' ~ 'ac/hc precursor',
                              CellType == 'Amacrine cells' ~ 'amacrine',
                              CellType == 'Bipolar cells' ~ 'bipolar',
                              CellType %in%  c('cones','Cones') ~ 'cone',
                              CellType == 'Corneal epithelium' ~ 'epithelium',
                              CellType %in% c('HCs','Horizontal cells') ~ 'horizontal',
                              CellType == 'Microglia' ~ 'microglia',
                              CellType == 'Muller glia' ~ 'mueller',
                              CellType == 'Photoreceptor precursors' ~ 'photoreceptor precursor',
                              CellType == 'RGCs' ~ 'retinal ganglion',
                              CellType %in% c("rods","Rods") ~ 'rod',
                              CellType %in%  c('RPC','RPCs') ~ 'rpc',
                              CellType == 'RPE' ~ 'rpe',
                              grepl("^T\\d", CellType) ~ 'neurogenic'
    )) %>% 
  mutate(SubCellType = CellType,
         CellType = MajorCellType,
         biosample_title = str_extract(file, "GSM\\d+"),
         bc = gsub('-\\d+', '', `...1`)) %>% 
  left_join(smeta_srp443999 %>% select(sample_accession, title, biosample_title), by = 'biosample_title') %>% 
  mutate(barcode = paste0(bc, '_', sample_accession)) %>% 
  select(barcode, MajorCellType, CellType, SubCellType, seurat_clusters, title)

sheets <- readxl::excel_sheets("data/41467_2024_47933_MOESM4_ESM.xlsx")
xcel_list <- list()
for (i in sheets[2:length(sheets)]){
  xcel_list[[i]] <- readxl::read_xlsx("data/41467_2024_47933_MOESM4_ESM.xlsx", sheet = i, skip = 3)
}

xcel_info <- xcel_list %>% bind_rows(.id = 'sheet') %>% select(sheet, seurat_clusters = cluster, author_excel_CT = `cell fate`) %>% 
  mutate(internal_id = gsub("PCW.*_","",sheet)) %>% 
  unique() %>% data.frame

srp443999_align <- srp443999_cells_cm %>% mutate(internal_id = gsub("_.*","",title)) %>% left_join(xcel_info, by = c("seurat_clusters", "internal_id")) 

# 15183_Retina_20PCW and 15184_Retina_21PCW are not correct. don't know why. just ignoring those two sheets.
srp443999_align %>% mutate(new_auth_CT = case_when(!title %in% c('15183_Retina_20PCW','15184_Retina_21PCW') & is.na(SubCellType) ~ author_excel_CT,
                                                   title %in% c('15183_Retina_20PCW','15184_Retina_21PCW') ~ SubCellType,
                                                   TRUE ~ SubCellType)) %>% 
  group_by(title, SubCellType, author_excel_CT, new_auth_CT) %>% summarise(Count = n()) %>% data.frame 

srp443999_align_corrected <- srp443999_align %>% mutate(new_auth_CT = case_when(!title %in% c('15183_Retina_20PCW','15184_Retina_21PCW') & is.na(SubCellType) ~ author_excel_CT,
                                                                                title %in% c('15183_Retina_20PCW','15184_Retina_21PCW') ~ SubCellType,
                                                                                TRUE ~ SubCellType)) %>% 
  mutate(SubCellType = new_auth_CT) %>%
  select(barcode, MajorCellType, CellType, SubCellType)

# Define a mapping of new labels to existing labels in ct.txt
label_map <- c(
  "Rods" = "rod",
  "Cones" = "cone",
  "RPCs" = "rpc",
  "Amacrine and horizontal cells" = NA,
  "HCs + ACs" = NA,
  "Mitotic cells" = NA,
  "Bipolar cells" = "bipolar",
  "RGCs" = "retinal ganglion",
  "periocular mesenchyme" = "periocular mesenchyme",
  "corneal stroma" = "cornea",
  "microglia" = "microglia",
  "corneal endothelium" = "endothelial",
  "melanocytes" = "melanocyte",
  "red blood cells" = "red blood",
  "ocular surface epithelium" = "epithelium",
  "NC" = "neural crest",
  "corneal nerves" = "corneal nerve",
  "trabecular mashwork" = "jct",
  "RPE" = "rpe",
  "extraocular muscles" = "muscle",
  "extraocular muscle" = "muscle",
  "periocular connective tissue" = "uveal",
  "Photoreceptor precursors" = "photoreceptor precursor",
  "Amacrine cells" = "amacrine",
  "Horizontal cells" = "horizontal",
  "Lens cells" = "fiber", 
  "Microglia" = "microglia",
  "Corneal epithelium" = "epithelium",
  "ACs" = "amacrine", 
  "rods" = "rod",
  "BPs" = "bipolar",
  "monocytes/macrophages" = "immune",
  "lens fibre" = "fiber",
  "corneal fibroblasts" = "fibroblast",
  "Muller glia" = "mueller",
  "RGC and amacrine cell precursor" = NA,
  "HCs" = "horizontal",
  "proliferating periocular connective tissue" = NA,
  "proliferating extraocular muscle" = "muscle",
  "proliferating cells" = NA,
  "RGC and amacrine cells" = NA,
  "Corneal stroma" = "cornea",
  "cones" = "cone",
  "blood vessel" = NA,
  "mixed cluster" = NA,
  "RPC" = "rpc",
  "Red blood cells" = "red blood",
  "corenal endothelium" = "endothelial",
  "cornea stroma" = "cornea",
  "corneal epithelium" = "epithelium"
)

# Apply relabeling
srp443999_align_corrected <- srp443999_align_corrected %>%
  mutate(cell_type = ifelse(SubCellType %in% names(label_map), label_map[SubCellType], SubCellType)) %>% 
  mutate(cell_type = case_when(grepl("^T\\d", cell_type) ~ "rpc",
                               TRUE ~ cell_type))

ct_out <- bind_rows(ct %>% filter(!barcode %in% srp443999_align_corrected$barcode),
          srp443999_align_corrected)
          
ct_out %>% write_csv('/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2025_02_28.01.csv.gz')
