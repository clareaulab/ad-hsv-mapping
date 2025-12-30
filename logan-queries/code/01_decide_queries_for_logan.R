library(data.table)
library(seqinr)
library(dplyr)
library(readxl)
library(BuenColors)

# Soh et al Nature Comms 2024 Table 2 = 41467_2024_54668_MOESM5_ESM.xlsx
# https://www.nature.com/articles/s41467-024-54668-2#Sec21

# Import Genbank data
hhv1 <- read.fasta("../../reference-genomes/HSV1-coding.fasta", as.string = TRUE, forceDNAtolower = FALSE )
hhv3 <- read.fasta("../../reference-genomes/VZV-coding.fasta", as.string = TRUE, forceDNAtolower = FALSE)

# Import data from nature comms paper
subset_2herpes <- readxl::read_excel("../data/41467_2024_54668_MOESM5_ESM.xlsx")[,c(4,5,7)] %>% data.frame()
colnames(subset_2herpes) <- c("what", "id_hhv1", "id_hhv3")
cc_subset_2herpes <- (subset_2herpes[complete.cases(subset_2herpes),])

# reprocess
hhv1_df <- lapply(hhv1, function(x){
  gene = stringr::str_split_fixed(attr(x, "Annot"), " ", 8)[,2]
  gene = gsub("]", "", gsub("[gene=", "", gene, fixed = TRUE),fixed = TRUE)
  data.frame(sequence = as.character(x), gene)
}) %>% rbindlist() %>% mutate(virus = "HHV1") %>% mutate(length = nchar(sequence)) %>%
  mutate(id_hhv1 = paste0("HSV-1_", gene))

hhv3_df <- lapply(hhv3, function(x){
  gene = stringr::str_split_fixed(attr(x, "Annot"), " ", 8)[,2]
  gene = gsub("]", "", gsub("[gene=", "", gene, fixed = TRUE),fixed = TRUE)
  data.frame(sequence = as.character(x), gene)
}) %>% rbindlist() %>% mutate(virus = "HHV3") %>% mutate(length = nchar(sequence)) %>%
  mutate(id_hhv3 = paste0("VZV_", gene))

# stich these all together
mdf <- merge(merge(cc_subset_2herpes, hhv1_df, by = "id_hhv1"), hhv3_df, by = "id_hhv3")

ggplot(mdf, aes(x = length.x, y = length.y)) + 
  geom_point() +
  scale_y_log10() + scale_x_log10() +
  labs(x = "Length of CDS - HSV-1", y = "Length of CDS - VZV") +
  pretty_plot(fontsize = 8) + L_border() 

(mdf %>% filter(length.x < 400 & length.y < 400))[,c("gene.x", "sequence.x")] %>% 
  write.table()

(mdf %>% filter(length.x < 400 & length.y < 400))[,c("gene.y", "sequence.y")] %>% 
  write.table()
