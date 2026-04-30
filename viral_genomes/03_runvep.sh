#!/bin/bash
#
#SBATCH --job-name=ranger
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --array=0-2

## get samples "ls /data1/lareauc/projects/logan_anecdotes/HHV1/fastqs/*R1_*   | cut -d / -f 8 | cut -d _ -f 1"
samples=( SRR19792154 SRR19792155 SRR19792156)

## make sure array is from 0:(n-1) samples
the_sample=${samples[$SLURM_ARRAY_TASK_ID]}
#the_sample=${samples[0]}

mkdir -p vep_outs

#http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer
vep_img="/data1/lareauc/users/gutierj/singularity_sifs/vep.sif"
the_path="/data1/lareauc/users/gutierj/projects_2024/gtex_hsv1/variant_calling"
singularity exec --bind .:/src $vep_img vep --dir $the_path -i /src/freebayes/${the_sample}.vcf -gff /src/hsv_ref/NC_001806.2.gff3.gz -fasta /src/hsv_ref/NC_001806.2.fa.gz -o /src/vep_outs/${the_sample}.vcf --force_overwrite

#singularity exec $vep_img ls ## look at all directories  
#singularity exec $vep_img cp ${the_sample}.vcf .
#singularity exec $vep_img cp ${the_sample}.vcf_summary.html .
#singularity exec $vep_img cp ${the_sample}.vcf_warnings.txt .




