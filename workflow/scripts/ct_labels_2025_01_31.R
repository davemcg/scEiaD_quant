library(tidyverse)
# add new studies

ct_labels <- data.table::fread('/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2024_08_27.csv.gz')


smeta <- data.table::fread('sample_meta.scEiaD_v1.2025_01_31.01.tsv.gz')
# pulling "temporary" full species h5ad objects made with ~/git/scEiaD_modeling/workflow/scripts/merge_adata.py
# to line up author samples <-> sra IDs by brute force (counting barcode matches)
mm_obs <- data.table::fread('data/mm111.adata.solo.2025x.obs.csv.gz') %>% select(-17)
hs_obs <- data.table::fread('~/data/scEiaD_modeling/hs111.adata.solo.2025x2.obs.csv.gz') %>% select(-17)


###################
# emtab9395
e9395_labels <- data.table::fread('data/E-MTAB-9395_Dev_path_width_Cell_stage_clique_umap_info.txt')

# multiple samples per time point - don't know what "1" and "2" refer to 
# have to align the barcodes to hand verify
mm_obs %>% filter(study_accession == 'E-MTAB-9395') %>% mutate(bc_core = gsub("_.*","",barcode)) %>% 
  left_join(e9395_labels %>% filter(source == 'Retinal') %>% mutate(bc_core = gsub('.*_','',V1)), by = 'bc_core') %>% 
  group_by(sample_accession, orig.ident) %>% summarise(Count = n()) %>% mutate(Ratio = Count/sum(Count)) %>% arrange(-Count) %>% 
  filter(!is.na(orig.ident))


e9395_labels_edits <- e9395_labels %>% 
  filter(source == 'Retinal') %>% 
  mutate(sample = case_when(orig.ident == 'Anouk_P5_1' ~ 'ERS4850920',
                            orig.ident == 'Anouk_P5_2' ~ 'ERS4850921',
                            orig.ident == 'Anouk_P9' ~ 'ERS4850922',
                            orig.ident == 'Anouk_P0_2' ~ 'ERS4850919',
                            orig.ident == 'Anouk_E13_1' ~ 'ERS4850916',
                            orig.ident == 'Anouk_E13_2' ~ 'ERS4850917',
                            orig.ident == 'Anouk_P0_1' ~ 'ERS4850918')) %>% 
  mutate(MajorCellType = case_when(grepl("rpc", tissue) ~ 'rpc',
                                   grepl("Precur", tissue) ~ 'photoreceptor precursor',
                                   tissue == 'NE' ~ 'neuroepithelium',
                                   tissue == 'RGCs' ~ 'retinal ganglion',
                                   tissue == 'Müller Glia' ~ 'mueller',
                                   tissue == 'RNeurogenic' ~ 'neurogenic',
                                   tissue == 'NI' ~ NA,
                                   TRUE ~ tissue)) %>% 
  mutate(MajorCellType = tolower(MajorCellType) %>% gsub("s$","",.) %>% gsub(" cell","",.),
         SubCellType = tissue) %>% 
  mutate(barcode = paste0(gsub('.*_','',V1), '_', sample))

mm_obs_9395 <- mm_obs %>% filter(study_accession == 'E-MTAB-9395') %>% 
  select(-MajorCellType, -SubCellType) %>% 
  left_join(e9395_labels_edits %>% select(barcode, MajorCellType, SubCellType), by = 'barcode')


#############
srp510712_cell <- data.table::fread('data/fetal_hrca_88444d73-7f55-4a62-bcfe-e929878c6c78.obs.csv.gz')


srp510712_relations <- hs_obs %>% filter(study_accession == 'SRP510712') %>% mutate(bc_core = gsub("_.*","",barcode) ) %>% 
  left_join(srp510712_cell %>% mutate(bc_core = gsub('.*_','',V1)%>% 
                                        gsub("-\\d","",.)), by = 'bc_core') %>% 
  mutate(sample = gsub("_[ATGC]{5,}.*","",V1)) %>% 
  group_by(sample_accession, sample) %>% summarise(Count = n()) %>% arrange(-Count) %>% 
  mutate(Ratio = Count/sum(Count)) %>% 
  filter(!is.na(sample)) %>% ungroup() %>% group_by(sample_accession) %>% 
  slice_max(n=1,order_by=Ratio) 

hs_obs_srp510712 <- hs_obs %>% filter(study_accession == 'SRP510712') %>% 
  mutate(bc_core = gsub("_.*","",barcode) ) %>% 
  select(-MajorCellType, -SubCellType, -cell_type_ontology_term_id) %>% 
  left_join(srp510712_relations, by = 'sample_accession') %>% 
  left_join(srp510712_cell %>% 
              mutate(bc_core = gsub('.*_','',V1)%>% 
                       gsub("-\\d","",.),
                     sample = gsub("_[ATGC]{5,}.*","",V1),
                     SubCellType = author_cell_type,
                     MajorCellType = case_when(grepl("Precursor", author_cell_type) ~ tolower(author_cell_type),
                                               grepl('RPC',author_cell_type) ~ 'rpc',
                                               grepl("bipolar|OFFx", cell_type) ~ 'bipolar',
                                               grepl("amacrine", cell_type) ~ "amacrine",
                                               grepl("horizontal", cell_type) ~ "horizontal",
                                               grepl("Mueller", cell_type) ~ "mueller",
                                               grepl("cone", cell_type) ~ "cone",
                                               grepl("retinal rod", cell_type) ~ "rod",
                                               grepl("ganglion", cell_type) ~ "retinal ganglion",
                     )), 
            by = c("bc_core","sample")) 


bind_rows(ct_labels,
          hs_obs_srp510712 %>% select(barcode, MajorCellType, CellType = MajorCellType, SubCellType, cell_type_ontology_term_id),
          mm_obs_9395 %>% select(barcode, MajorCellType, CellType = MajorCellType, SubCellType)) %>% 
  write_csv(, file = '/Users/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2025_01_31.01.csv.gz')
