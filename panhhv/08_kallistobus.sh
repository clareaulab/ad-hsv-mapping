

#conda activate /data1/lareauc/users/gutierj/env/kb


# Build the index and generate the t2g file
## This is kallisto: 0.52.0 index 
#kb ref \
#    -i hg38_panhhv_kallisto/hg38_panhhv_trans.idx \
#    -g hg38_panhhv_kallisto/hg38_panhhv_t2g.txt \
#    -f1 hg38_panhhv_kallisto/hg38_panhhv_transcriptome.fa \
#    hg38_panhhv.fasta \
#    hg38_panhhv.gtf

mkdir -p hg38_panhhv_kallisto_hardmasked
kb ref \
   -i hg38_panhhv_kallisto_hardmasked/hg38_panhhv_trans.idx \
      -g hg38_panhhv_kallisto_hardmasked/hg38_panhhv_t2g.txt \
    -f1 hg38_panhhv_kallisto_hardmasked/hg38_panhhv_transcriptome.fa \
    hardmasked_hg38_panhhv.fasta \
    hg38_panhhv.gtf



