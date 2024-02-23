#!/bin/bash
#SBATCH --partition normal
#SBATCH --job-name=starRNA
#SBATCH --output=starRNA_%j
#SBATCH --mem=50G
#SBATCH -n 24
#SBATCH --time=48:0:0

#maybe run ~/.bashrc?

source ~/.bashrc

conda activate bio-python2
module load STAR/2.7.11a 
module load samtools/1.9

REF=/gs/gsfs0/users/kpradha1/projects/star_hg38


SNAME=$1
FQ1=$2
FQ2=$3


ulimit -n 10000


STAR --genomeDir $REF \
--runThreadN 24 \
--readFilesIn $FQ1 $FQ2 \
--outFileNamePrefix $SNAME/ \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType                  BySJout \
--outFilterMultimapNmax          20      \
--alignSJoverhangMin             8       \
--alignSJDBoverhangMin           1       \
--outFilterMismatchNmax          999     \
--outFilterMismatchNoverLmax     0.06    \
--alignIntronMin                 20      \
--alignIntronMax                 1000000 \
--alignMatesGapMax               1000000 \
--twopassMode None \
--quantMode GeneCounts

#rename file and index it
cp $SNAME/Aligned.sortedByCoord.out.bam $SNAME/${SNAME}_sorted.bam
samtools index $SNAME/${SNAME}_sorted.bam $SNAME/${SNAME}_sorted.bai

#gene counts
cp $SNAME/ReadsPerGene.out.tab $SNAME/${SNAME}_ReadsPerGene.out.tab
