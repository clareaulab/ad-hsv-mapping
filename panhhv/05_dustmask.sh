#!/bin/bash
#
#SBATCH --job-name=cluster_nr
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G



module load lareauc/blast/2.16
module load lareauc/bedtools/2.30

in_fa="../pan_hhv_genomic.fna" 

## Move into directory to make cleaner 
mkdir -p dustmask
cd dustmask

dustmasker -window 30 -in ${in_fa} -outfmt interval > dust30.int

# Convert interval dust file to bed file
while read -r line; do 
	if [ "$(echo "$line" | head -c 1)" = ">" ]; then                                 
	#header=$(echo "$line" | sed 's/>//g' - )
	header=$(echo "$line" | sed -n 's/.*>\([^ ]*\).*/\1/p' - ) 
	else
	start=$(echo "$line" | cut -f1 -d' ' -)
	end=$(echo "$line" | cut -f3 -d' ' - )
	echo -e "${header}\t${start}\t${end}" >> dust30.bed
	fi 
done < dust30.int



## Now do 64 
dustmasker -window 64 -in ${in_fa} -outfmt interval > dust64.int

while read -r line; do 
	if [ "$(echo "$line" | head -c 1)" = ">" ]; then                                 
	#header=$(echo "$line" | sed 's/>//g' - )
	header=$(echo "$line" | sed -n 's/.*>\([^ ]*\).*/\1/p' - ) 
	else
	start=$(echo "$line" | cut -f1 -d' ' -)
	end=$(echo "$line" | cut -f3 -d' ' - )
	echo -e "${header}\t${start}\t${end}" >> dust64.bed
	fi 
done < dust64.int



# Merge the bed files
cat dust30.bed dust64.bed > mask.tmp


# Sort the cat bed file
sort -k1,1 -k2,2n mask.tmp > mask.sort.tmp

# Clean up some bugs (-1 and a space)
## IDK WHY??? 
sed 's/ //g' mask.sort.tmp \
	  | sed 's/-1/0/g' - \
	    > mask.sort.clean.tmp


bedtools merge -i mask.sort.clean.tmp > mask.regions.tmp


# Clean up some bugs (-1 and a space)
sed 's/ //g' mask.regions.tmp \
	  | sed 's/-1/0/g' - \
	    > mask.regions.bed

bedtools maskfasta -fi $in_fa -bed mask.regions.bed -fo hardmasked_pan_hhv_genomic.fna

bedtools maskfasta -fi $in_fa -bed mask.regions.bed -fo softmasked_pan_hhv_genomic.fna -soft


cp hardmasked_pan_hhv_genomic.fna ../
