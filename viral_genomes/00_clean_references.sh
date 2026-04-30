

conda activate /data1/lareauc/users/gutierj/env/htslib
#ml samtools

mkdir -p hsv_ref

bgzip -c ../NC_001806.2.fasta > hsv_ref/NC_001806.2.fa.gz

#tabix  hsv_ref/NC_001806.2.fa.gz


##http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff
grep -v "#" ../NC_001806.2.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > hsv_ref/NC_001806.2.gff3.gz
tabix -p gff hsv_ref/NC_001806.2.gff3.gz


