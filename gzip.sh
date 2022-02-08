#!/bin/bash
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuozhen.bao@yale.edu
#SBATCH --job-name=rezip
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=16g 
#SBATCH --time=120:00:00

gzip -d hC2_CKDL210027950-1a_H2YCKDSX3_L1_2.fq.gz
gzip -d hC2_CKDL210027950-1a_H2YCKDSX3_L1_1.fq.gz

gzip hC2_CKDL210027950-1a_H2YCKDSX3_L1_2.fq
gzip hC2_CKDL210027950-1a_H2YCKDSX3_L1_1.fq

