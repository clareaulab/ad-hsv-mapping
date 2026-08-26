#conda activate /data1/lareauc/users/gutierj/env/entrez

readarray -t the_accs < pan_hhv_accs.txt


mkdir -p gffs/


for acc in ${the_accs[@]}; do 
#acc=NC_001806.2
efetch -db nuccore -id "${acc}" -format gff3 > gffs/${acc}.gff3

done
