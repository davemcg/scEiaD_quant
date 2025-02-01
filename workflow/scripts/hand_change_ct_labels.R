library(tidyverse)
cmeta <- data.table::fread("~/git/scEiaD_quant/scEiaD_cell_labels_2024_03_07.csv.gz")

cmeta2 <- cmeta %>% mutate(MajorCellType = case_when(grepl("margin/peri", MajorCellType) ~ NA,
                                           MajorCellType == 'acta2' ~ 'smooth muscle',
                                           MajorCellType %in% c('neutrophil','smc','tuft') ~ NA,
                                           MajorCellType == 'doublet' ~ NA,
                                           grepl("endothelial", MajorCellType) ~ "endothelial",
                                           grepl("epithelial", MajorCellType) ~ "epithelial",
                                           MajorCellType == 'earlyfiber' ~ 'fiber',
                                           MajorCellType == 'myeschwann' ~ 'schwann',
                                           grepl("rpc", MajorCellType) ~ 'rpc',
                                           MajorCellType == 'blood vessel' ~ NA,
                                           MajorCellType == 'transitional' ~ 'epithelial',
                                           MajorCellType == 'ciliary muscle' ~ 'muscle',
                                           TRUE ~ MajorCellType
                                           ))

cmeta2 %>% write_csv("~/git/scEiaD_quant/scEiaD_cell_labels_2024_08_27.csv.gz")


