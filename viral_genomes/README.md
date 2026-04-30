Caleb asked me to do variant effect prediction on the HSV1 data. 

https://github.com/freebayes/freebayes

0) samtools view HSV1 only 
1) freebayes on the HSV1 reads using HSV1 reference.fa
2) VEP using the VCFs from ^^^ and the .gtf/.gff3 filei

https://www.ensembl.info/2020/08/28/cool-stuff-the-ensembl-vep-can-do-annotating-sars-cov-2-variants/
https://www.ensembl.info/2018/07/26/cool-stuff-the-vep-can-do-custom-annotation/


VEP install is painfull I just am using the singularity image... hopefully I can just deploy it easily 
http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer

Here: `/data1/lareauc/users/gutierj/singularity_sifs/vep.sif`

freebayes -f ref.fa aln.bam >var.vcf

./vep -i input.vcf -gff data.gff.gz -fasta genome.fa.gz 


