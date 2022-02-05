#!/bin/bash
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuozhen.bao@yale.edu
#SBATCH --job-name=st-pipeline
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem-per-cpu=4g 
#SBATCH --time=120:00:00

sh stpipeline.sh hC2 /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/02.effective
