#!/bin/bash
#SBATCH --job-name=Admixture
#SBATCH --time=30-00:00:00
#SBATCH --partition=unlimited
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --output=ADMIX.%A_%a.out

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate omics-workshop-python 

OUT="admix_results"

mkdir -p $OUT

for K in {2..20}; do 
	admixture --cv admixture_panel.thin.bed  $K | tee $OUT/log$K.out; 
done
