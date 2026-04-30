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

bayes=/data1/lareauc/users/gutierj/bin/freebayes-1.3.6-linux-amd64-static

mkdir -p freebayes

the_ref="../NC_001806.2.fasta"

$bayes -f $the_ref viral_bams/${the_sample}_hsv1.bam > freebayes/${the_sample}.vcf


