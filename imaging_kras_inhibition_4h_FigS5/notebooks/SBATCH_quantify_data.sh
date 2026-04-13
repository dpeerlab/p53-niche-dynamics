#!/bin/bash
#SBATCH --job-name=quantify
#SBATCH --time=3:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mem=128G
conda activate spatialdata_opencv
python Krasi4h_03_quantify_channels_batch.py --slide JR-sp-15-010
