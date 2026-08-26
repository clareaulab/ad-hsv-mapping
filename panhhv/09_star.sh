
gtf="hg38_panhhv.gtf"
fasta="hg38_panhhv.fasta"
# Generate STAR reference
/usersoftware/lareauc/STAR_2.7.11b/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN 32 --genomeDir hg38_panhhv_STAR_reference --genomeFastaFiles $fasta --sjdbGTFfile $gtf

