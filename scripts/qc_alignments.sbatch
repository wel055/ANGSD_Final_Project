#!/bin/bash
#SBATCH --job-name=QC_alignments
#SBATCH --output=logs/qc_align_%j.out
#SBATCH --error=logs/qc_align_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=angsd_class

echo "Running on compute node: $(hostname)"
echo "Start time: $(date)"

source /home/fs01/wel4014/miniforge3/etc/profile.d/conda.sh
conda activate angsd

PROJECT_DIR=/home/fs01/wel4014/angsd_proj
ALIGN_DIR=${PROJECT_DIR}/alignments
FASTQC_DIR=${PROJECT_DIR}/qc/fastqc
MULTIQC_DIR=${PROJECT_DIR}/qc/multiqc

mkdir -p ${FASTQC_DIR}
mkdir -p ${MULTIQC_DIR}

echo "===================================="
echo "Running FastQC on alignment BAM files"
echo "===================================="

fastqc \
    --outdir ${FASTQC_DIR} \
    --threads 4 \
    ${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam

echo "FastQC completed: $(date)"

echo "===================================="
echo "Running MultiQC to aggregate results"
echo "===================================="

conda deactivate
conda activate multiqc

multiqc \
    --outdir ${MULTIQC_DIR} \
    --filename qc_alignments_report \
    --force \
    ${FASTQC_DIR}

echo "MultiQC completed: $(date)"
echo "Done! Check ${MULTIQC_DIR}/qc_alignments_report.html"
