library(dplyr)

args = commandArgs(trailingOnly=TRUE)

geneids <- read.csv(args[1], header = F)
genenames <- read.csv(args[2], header = F)

gene_tab <- cbind(geneids, genenames)
colnames(gene_tab) <-  c('id','name')

gene_tab <- gene_tab  %>% mutate(nname = case_when(id != name ~ paste0(name, ' (', id, ')'), TRUE ~ id))

write.csv(gene_tab, file = args[3])

