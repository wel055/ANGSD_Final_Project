#!/bin/bash
#SBATCH --job-name=STAR_4samples
#SBATCH --output=logs/star_align4_%j.out
#SBATCH --error=logs/star_align4_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=angsd_class

echo "Running on compute node: $(hostname)"
echo "Start time: $(date)"

source /home/fs01/wel4014/miniforge3/etc/profile.d/conda.sh
conda activate angsd

GENOME_DIR=/home/fs01/wel4014/angsd_proj/star_index
PROJECT_DIR=/home/fs01/wel4014/angsd_proj
ALIGN_DIR=${PROJECT_DIR}/alignments

mkdir -p ${ALIGN_DIR}

# New control samples
CONTROL_SAMPLES=("SRR5305483" "SRR5305576")

# New AD samples
AD_SAMPLES=("SRR5305558" "SRR5305477")

ALL_SAMPLES=("${CONTROL_SAMPLES[@]}" "${AD_SAMPLES[@]}")

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    echo "===================================="
    echo "Aligning ${SAMPLE}"

    # Determine directory
    if [[ " ${CONTROL_SAMPLES[@]} " =~ " ${SAMPLE} " ]]; then
        FASTQ_DIR=${PROJECT_DIR}/raw_sra/control
    else
        FASTQ_DIR=${PROJECT_DIR}/raw_sra/AD
    fi

    R1=${FASTQ_DIR}/${SAMPLE}_1.fastq.gz
    R2=${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

    STAR \
        --runThreadN 8 \
        --genomeDir ${GENOME_DIR} \
        --readFilesIn ${R1} ${R2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${ALIGN_DIR}/${SAMPLE}. \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic

    samtools index ${ALIGN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam

    samtools flagstat ${ALIGN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam \
        > ${ALIGN_DIR}/${SAMPLE}.flagstat.txt

done

echo "Alignment completed: $(date)"
