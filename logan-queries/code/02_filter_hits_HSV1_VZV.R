library(data.table)
library(dplyr)
library(BuenColors)

filter_logan_query <- function(df){
  df %>% filter(organism == "Homo sapiens") %>% 
    filter(!grepl("HSV|viral", bioproject_description)) %>% # remove known infections
    filter(!grepl("rganoid|PSC", bioproject_description)) %>% # remove culture models
    filter(grepl("ementia",bioproject_description)) # require brain/ dementia
}

# HSV-1 query
(fread("../logan-search-output/HHV1_UL31.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/HHV1_UL11.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/HHV1_US9.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/HHV1_UL49A.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]

# VZV query
(fread("../logan-search-output/VZV_ORF49.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/VZV_ORF65.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/VZV_ORF9A.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
