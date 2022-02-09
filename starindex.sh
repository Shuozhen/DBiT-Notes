#!/bin/bash
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuozhen.bao@yale.edu
#SBATCH --job-name=genomeGenerate
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem-per-cpu=4g 
#SBATCH --time=120:00:00

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /gpfs/ysm/project/fan/sb2723/00.database/hg38/STARindex --genomeFastaFiles /gpfs/ysm/project/fan/sb2723/00.database/hg38/hg38.fa --sjdbGTFfile /gpfs/ysm/project/fan/sb2723/00.database/hg38/gencode.v39.annotation.gtf --sjdbOverhang 149
