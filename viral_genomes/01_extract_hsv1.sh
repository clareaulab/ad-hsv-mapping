#!/bin/bash
#
#SBATCH --job-name=subset
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

the_samples=( SRR19792154 SRR19792155 SRR19792156 )

ml samtools 

mkdir -p viral_bams

parallel 'samtools view -b ../{}/outs/possorted_genome_bam.bam NC_001806.2 > viral_bams/{}_hsv1.bam' ::: ${the_samples[@]}


