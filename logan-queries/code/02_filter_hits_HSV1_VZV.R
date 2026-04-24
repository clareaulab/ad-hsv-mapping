library(data.table)
library(dplyr)
library(BuenColors)

filter_logan_query <- function(df){
  df %>% filter(organism == "Homo sapiens") %>% 
    filter(!grepl("HSV|viral|nfect", bioproject_description)) %>% # remove known infections
    filter(!grepl("rganoid|PSC|mouse", bioproject_description)) %>% # remove culture models
    filter(grepl("ementia|Brain|brain|lzhei",bioproject_description)) # require brain/ dementia
}

# HSV-1 query
(fread("../logan-search-output/HHV1_UL11.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/HHV1_US9.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/HHV1_UL49A.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]

ids_UL11 <- fread("../logan-search-output/HHV1_UL11.tsv") %>% filter(organism == "Homo sapiens") %>% pull(ID)
ids_US9 <- fread("../logan-search-output/HHV1_US9.tsv") %>% filter(organism == "Homo sapiens") %>% pull(ID)
ids_L49A<- fread("../logan-search-output/HHV1_UL49A.tsv") %>% filter(organism == "Homo sapiens") %>% pull(ID)

get_venn_counts_df <- function(vecA, vecB, vecC) {
  # Convert to unique sets
  A <- unique(vecA)
  B <- unique(vecB)
  C <- unique(vecC)
  
  # 1. Triple Intersection (The Center)
  abc <- intersect(intersect(A, B), C)
  
  # 2. Double Intersections (Excluding the triple center)
  ab_only <- setdiff(intersect(A, B), abc)
  bc_only <- setdiff(intersect(B, C), abc)
  ac_only <- setdiff(intersect(A, C), abc)
  
  # 3. Exclusive Segments (Unique to one set only)
  a_only <- setdiff(A, union(B, C))
  b_only <- setdiff(B, union(A, C))
  c_only <- setdiff(C, union(A, B))
  
  # Create the final dataframe
  data.frame(
    Region = c("A_only", "B_only", "C_only", "AB_only", "BC_only", "AC_only", "ABC"),
    Count = c(length(a_only), length(b_only), length(c_only), 
              length(ab_only), length(bc_only), length(ac_only), 
              length(abc)),
    stringsAsFactors = FALSE
  )
}
get_venn_counts_df(ids_UL11,ids_US9,ids_L49A)

get_venn_counts_df(
   fread("../logan-search-output/HHV1_UL11.tsv") %>% filter_logan_query()  %>% pull(ID),
  fread("../logan-search-output/HHV1_US9.tsv") %>% filter_logan_query()  %>% pull(ID),
  fread("../logan-search-output/HHV1_UL49A.tsv") %>% filter_logan_query() %>% pull(ID)
  
)


# VZV query
(fread("../logan-search-output/VZV_ORF49.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/VZV_ORF65.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]
(fread("../logan-search-output/VZV_ORF9A.tsv") %>% filter_logan_query)[,c("ID", "kmer_coverage", "bioproject_title")]


