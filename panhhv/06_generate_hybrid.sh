#!/bin/bash
#
#SBATCH --job-name=mkref
#SBATCH --partition=cpu,lareauc_cpu
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

## Cell Ranger hg38

HUMAN_FASTA_PATH="/data1/lareauc/references/txg/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
VIRAL_FASTA_PATH="hardmasked_pan_hhv_genomic.fna"

HUMAN_GTF_PATH="/data1/lareauc/references/txg/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
VIRAL_GTF_PATH="clean_pan_hhv_combined.gtf"


#COMBINED_FASTA_PATH="hg38_panhhv.fasta"
COMBINED_FASTA_PATH="hardmasked_hg38_panhhv.fasta"

COMBINED_GTF_PATH="hg38_panhhv.gtf"

cat ${HUMAN_FASTA_PATH} ${VIRAL_FASTA_PATH} > ${COMBINED_FASTA_PATH}


#####cat ${HUMAN_GTF_PATH} ${VIRAL_GTF_PATH} > ${COMBINED_GTF_PATH}
zcat ${HUMAN_GTF_PATH}  > ${COMBINED_GTF_PATH}
cat ${VIRAL_GTF_PATH} >>  ${COMBINED_GTF_PATH}

