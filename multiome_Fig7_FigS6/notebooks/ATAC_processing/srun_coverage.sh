#!/bin/bash

# Use ATAC_processing directory as input
f=$(pwd)
proj=$f/multiome_PP20251110/MotifData/MotifInstances
for i in `ls $proj`
do
filename_without_ext=$(basename "$i" ".tsv.gz")
srun --job-name=$filename_without_ext --time=3:00:00 --mem=16G --partition=cpu Rscript --quiet --vanilla --no-save $f/step09_compute_motif_coverage.R multiome_PP20251110 SEACell $proj/$i 1 147 > $filename_without_ext.stdout &
done
