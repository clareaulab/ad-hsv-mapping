
## RUN THiS ON INTERACTIVE ESSESION
## with CPU only took more than 20 min. 
## with 4 GPU cores it was <1 min SLAY!
alias isgpustrunk="salloc -p lareauc_gpu,gpu,gpushort --gres=gpu:4 --mem 32G -t 2:00:00 -c 16 -J isgpustrunk"

esm_path="/home/gutierj6/gutierj/bin/esm-variants"

python ${esm_path}/esm_score_missense_mutations.py --input-fasta-file NC_001806.2.faa --output-csv-file NC_001806.2_esm_scores.csv

