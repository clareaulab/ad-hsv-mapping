
# NCBI Command Line Tool 
#conda activate /data1/lareauc/users/gutierj/env/entrez

#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=node&id=3044472&lvl=
use_tax=3044472 ## Orthoherpesviridae/Herpesviridae

## Round 1: Download metadata for relevant references
datasets summary virus genome taxon ${use_tax} --complete-only --host human --refseq --as-json-lines | dataformat tsv virus-genome > herpesviridae_summary.tsv


## Now download Fasta and proteins
#datasets download virus genome taxon ${use_tax} --complete-only --host human --include genome,cds,protein,annotation --filename herpesviridae.zip 

## Extracting the 10 hhv 
awk -F'\t' 'NR>1 {print $1}'  herpesviridae_summary.tsv > raw_pan_hhv_accs.txt

## Of course HSV1 is not in this query as it is not annotated as host==human 
echo "NC_001806.2"  >> raw_pan_hhv_accs.txt


## Manually curate final genomes using vclust, exclude any pairs with >95 ANI and take longest genome. 
## 2 EBV genomes 98% ANI taking NC_007605.1 due to better annotations
##  NC_009334.1 is longer but unannotated
## HHV6A and HHV6B are 90% ANI so I will keep them seperate 
## HHV1 and HHV2 are 70%  ANI
cat raw_pan_hhv_accs.txt | grep -v NC_009334.1 > pan_hhv_accs.txt

## Download the summary info again
datasets summary virus genome accession --inputfile pan_hhv_accs.txt --as-json-lines | dataformat tsv virus-genome > human_herpesviridae_summary.tsv

## Download the fasta info
datasets download virus genome accession --inputfile pan_hhv_accs.txt --include genome,cds,protein,annotation --filename human_herpesviridae.zip



