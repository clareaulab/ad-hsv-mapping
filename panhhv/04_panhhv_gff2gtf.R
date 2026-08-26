#conda activate /data1/lareauc/users/gutierj/env/rbio
library(rtracklayer);library(dplyr) 
check_genecnt <- function(gff3) gff3 %>% as.data.frame()  %>% count(Name)  %>% filter(n>1) 


#dedup_gene_ids <- function(gr) {
#	  gid <- dplyr::coalesce(gr$gene_id, gr$gene, gr$Name, gr$ID, "unassigned_feature")
#  runs <- rle(as.character(gid))
#    runs$values <- make.unique(runs$values)
#    gr$gene_id <- rep(runs$values, runs$lengths)
#      gr$transcript_id <- gr$gene_id
#      return(gr)
#}

dedup_gene_ids <- function(gr) {
  gid <- dplyr::coalesce(gr$gene_id, gr$gene, gr$Name, gr$ID, "unassigned_feature")
  gid <- paste(as.character(seqnames(gr)), gid, sep="_")
  runs <- rle(gid)
  runs$values <- make.unique(runs$values)
  gr$gene_id <- rep(runs$values, runs$lengths)
  gr$transcript_id <- gr$gene_id
  return(gr)
}

gff_l <- list() ## final list I want to export

## HSV1 
aref <- "NC_001806.2"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## HSV1 has an inverted repeat region that is found in 2 genomic locations so deleting one and so all reads map to the second
## In fact there are two terminal repeats, masking the ends to everything piles up in the middle
gff3 <- gff3[( end(gff3) > 9213 & end(gff3) < 146737 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))

export(gff3,agtf,"gtf")



## HSV2
aref <- "NC_001798.2"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## HSV2 also has the same structure different coordinates obvi
gff3 <- gff3[( end(gff3) > 9300 & end(gff3) < 148034 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")



## VSV
aref <- "NC_001348.1"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## VZV has a single terminal repeat on the 3' end 
gff3 <- gff3[( end(gff3) < 117585 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")


## EBV 
aref <- "NC_007605.1"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)
span <- end(gff3) > 171823 & start(gff3) <= 171823
over <- start(gff3) > 171823
dup <- gff3[span]
start(dup) <- 1
end(dup) <- end(dup) - 171823
end(gff3)[span] <- 171823
start(gff3)[over] <- start(gff3)[over] - 171823
end(gff3)[over] <- end(gff3)[over] - 171823
gff3 <- c(gff3, dup)

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l, list(gff3))
export(gff3, agtf, "gtf")


## CMV
aref <- "NC_006273.2"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## No genes in repeat regions but duplciate genes potentially double check 
gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")


## HHV6A
aref <- "NC_001664.4"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## DRS repeats on 5' and 3' end, convention says to mask the 3' end 
gff3 <- gff3[( end(gff3) < 151235 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")



## HHV6B
aref <- "NC_000898.1"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## Again, DRS both sides masking 3'
gff3 <- gff3[( end(gff3) < 153322 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")


## HHV7
aref <- "NC_001716.2"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## Same as before masking 3' end 
gff3 <- gff3[( end(gff3) < 143047 )]

gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")


## HHV8
aref <- "NC_009333.1"
agff <- sprintf("gffs/%s.gff3",aref)
agtf <- sprintf("gffs/%s.gtf",aref)
gff3 <- import(agff)

## NO repeats in genome clear! 
gff3 <- dedup_gene_ids(gff3)

gff_l <- c(gff_l,list(gff3))
export(gff3,agtf,"gtf")



## COMBINE EVERYTHING
gff_comb <- do.call(c,gff_l)
export(gff_comb, "pan_hhv_combined.gtf", format="gtf")


gff_comb$type[gff_comb$type == "CDS"] <- "exon"
has_exon <- gff_comb$transcript_id[gff_comb$type == "exon"]
orphan_mask <- !(gff_comb$transcript_id %in% has_exon) & !(gff_comb$type %in% c("gene", "region"))
gff_comb$type[orphan_mask] <- "exon"
#gff_comb$gene_id <- dplyr::coalesce(gff_comb$gene, gff_comb$Name, gff_comb$ID, "unassigned_feature")
#gff_comb$transcript_id <- gff_comb$gene_id


df <- as.data.frame(gff_comb) %>% 
	  dplyr::group_by(seqnames, gene_id) %>% 
	    dplyr::summarize(start = min(start), end = max(end), .groups = "drop")
    write.table(df, "viral_gene_coordinates.tsv", sep="\t", quote=FALSE, row.names=FALSE)

rtracklayer::export(gff_comb, "clean_pan_hhv_combined.gtf", format="gtf")

