#!/bin/bash
#SBATCH --partition normal
#SBATCH --job-name=peakPipeline
#SBATCH --output=peakPipeline_%j
#SBATCH --mem=50G
#SBATCH -n 1
#SBATCH --ntasks-per-node=10
#SBATCH --time=48:0:0

#maybe run ~/.bashrc?

source ~/.bashrc

conda activate bio-python2

PATH=$PATH:/gs/gsfs0/hpc01/rhel8/apps/conda3/envs/biobase/bin/

REF=~/projects/GRCh38_no_alt_analysis/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index

SNAME=$1
FQ1=$2
FQ2=$3


bowtie2 -t -p 8 -q --local -x $REF -1 <(zcat $FQ1) -2 <(zcat $FQ2) 2>| \
    summary_${SNAME}.txt | samtools view -bS - >| ${SNAME}.bam

picard SortSam INPUT=${SNAME}.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate \
    CREATE_INDEX=true MAX_RECORDS_IN_RAM=1000000 TMP_DIR=tmp

picard MarkDuplicates INPUT=${SNAME}_sorted.bam OUTPUT=${SNAME}_sorted_mdup.bam \
    METRICS_FILE=${SNAME}_mdMetric.txt 

picard BuildBamIndex INPUT=${SNAME}_sorted_mdup.bam



BAM=${SNAME}_sorted_mdup.bam
FOL=macs_${SNAME}
macs2 callpeak -f BAMPE -t $BAM -n $SNAME --outdir $FOL --verbose 3



