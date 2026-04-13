#!/bin/bash
#SBATCH --job-name=quantify
#SBATCH --time=10:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mem=128G
conda activate spatialdata_opencv
python ProgenitorNiche_03_quantify_channels_batch.py --slide SLIDE-0289
