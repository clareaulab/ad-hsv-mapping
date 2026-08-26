
## Next we extract the data retaining genomic,cds, and proteins
unzip human_herpesviridae.zip

cp ncbi_dataset/data/genomic.fna pan_hhv_genomic.fna
cp ncbi_dataset/data/cds.fna pan_hhv_cds.fna
cp ncbi_dataset/data/protein.faa pan_hhv_protein.faa

rm -rf ncbi_dataset/


