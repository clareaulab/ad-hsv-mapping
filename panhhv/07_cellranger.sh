#!/bin/bash
#
#SBATCH --job-name=mkref
#SBATCH --partition=cpu,lareauc_cpu
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

## for iris
 ml cellranger

#/data1/lareauc/references/txg/refdata-gex-GRCm39-2024-A/refdata-gex-GRCh38-2024-A


COMBINED_FASTA_PATH="hardmasked_hg38_panhhv.fasta"

COMBINED_GTF_PATH="hg38_panhhv.gtf"

## Test quick alignment to confirm its working 
#cellranger mkref --genome=hsv1_test --fasta=${VIRAL_FASTA_PATH} --genes=${VIRAL_GTF_PATH}

cellranger mkref --genome=hg38_panhhv_hardmasked_cellranger --nthreads=32 --memgb=64 \
    --fasta=${COMBINED_FASTA_PATH} \
    --genes=${COMBINED_GTF_PATH}

