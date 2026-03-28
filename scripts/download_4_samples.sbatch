#!/bin/bash
#SBATCH --job-name=SRA_download
#SBATCH --output=logs/sra_download_%j.out
#SBATCH --error=logs/sra_download_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=angsd_class

echo "Running on compute node: $(hostname)"
echo "Start time: $(date)"

source /home/fs01/wel4014/miniforge3/etc/profile.d/conda.sh
conda activate angsd

PROJECT_DIR=/home/fs01/wel4014/angsd_proj

# AD samples
AD_SAMPLES=("SRR5305558" "SRR5305477")
# Control samples
CONTROL_SAMPLES=("SRR5305483" "SRR5305576")

# Download AD samples
for SAMPLE in "${AD_SAMPLES[@]}"; do
    echo "===================================="
    echo "Downloading ${SAMPLE} (AD)"
    echo "===================================="
    fasterq-dump --split-files --threads 4 \
        --outdir ${PROJECT_DIR}/raw_sra/AD \
        ${SAMPLE}
    echo "Compressing ${SAMPLE}..."
    gzip ${PROJECT_DIR}/raw_sra/AD/${SAMPLE}_1.fastq
    gzip ${PROJECT_DIR}/raw_sra/AD/${SAMPLE}_2.fastq
    echo "${SAMPLE} done: $(date)"
done

# Download Control samples
for SAMPLE in "${CONTROL_SAMPLES[@]}"; do
    echo "===================================="
    echo "Downloading ${SAMPLE} (Control)"
    echo "===================================="
    fasterq-dump --split-files --threads 4 \
        --outdir ${PROJECT_DIR}/raw_sra/control \
        ${SAMPLE}
    echo "Compressing ${SAMPLE}..."
    gzip ${PROJECT_DIR}/raw_sra/control/${SAMPLE}_1.fastq
    gzip ${PROJECT_DIR}/raw_sra/control/${SAMPLE}_2.fastq
    echo "${SAMPLE} done: $(date)"
done

echo "All downloads completed: $(date)"
