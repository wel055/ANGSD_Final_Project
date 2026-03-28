#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=logs/featurecounts_%j.out
#SBATCH --error=logs/featurecounts_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=angsd_class

echo "Running on compute node: $(hostname)"
echo "Start time: $(date)"

source /home/fs01/wel4014/miniforge3/etc/profile.d/conda.sh
conda activate angsd

PROJECT_DIR=/home/fs01/wel4014/angsd_proj
ALIGN_DIR=${PROJECT_DIR}/alignments
GTF=/athena/angsd/scratch/yiw4021/project/refseq/GCF_000001405.40_GRCh38.p14_genomic.gtf
OUT_DIR=${PROJECT_DIR}/counts

mkdir -p ${OUT_DIR}

# All 10 samples
BAMS="${ALIGN_DIR}/SRR5305524.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305480.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305508.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305546.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305567.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305591.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305483.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305576.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305558.Aligned.sortedByCoord.out.bam \
${ALIGN_DIR}/SRR5305477.Aligned.sortedByCoord.out.bam"

echo "===================================="
echo "Running featureCounts (gene-level)"
echo "===================================="

featureCounts \
    -a ${GTF} \
    -o ${OUT_DIR}/gene_counts.txt \
    -g gene_id \
    -t exon \
    -T 8 \
    -p --countReadPairs \
    -s 0 \
    ${BAMS}

echo "featureCounts completed: $(date)"
