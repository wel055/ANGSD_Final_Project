#!/bin/bash
#SBATCH --job-name=QC_align4
#SBATCH --output=logs/qc_align4_%j.out
#SBATCH --error=logs/qc_align4_%j.err
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

# Run FastQC only on the 4 new BAMs
NEW_SAMPLES=("SRR5305558" "SRR5305477" "SRR5305483" "SRR5305576")

echo "===================================="
echo "Running FastQC on new alignment BAM files"
echo "===================================="

for SAMPLE in "${NEW_SAMPLES[@]}"; do
    fastqc \
        --outdir ${FASTQC_DIR} \
        --threads 4 \
        ${ALIGN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam
done

echo "FastQC completed: $(date)"

echo "===================================="
echo "Running MultiQC on ALL samples (original 6 + new 4)"
echo "===================================="

conda deactivate
conda activate multiqc

# Run MultiQC on entire fastqc dir (includes all 10 samples)
multiqc \
    --outdir ${MULTIQC_DIR} \
    --filename qc_all_10_samples \
    --force \
    ${FASTQC_DIR}

echo "MultiQC completed: $(date)"
echo "Done! Check ${MULTIQC_DIR}/qc_all_10_samples.html"
