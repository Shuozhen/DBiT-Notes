#!/bin/bash
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuozhen.bao@yale.edu
#SBATCH --job-name=genomeGenerate_nc
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem-per-cpu=4g 
#SBATCH --time=120:00:00


STAR --runThreadN 16 --limitGenomeGenerateRAM 50000000000 --genomeSAindexNbases 12 --runMode genomeGenerate --genomeDir /gpfs/ysm/project/fan/sb2723/00.database/hg38/StarIndex_nc --genomeFastaFiles /gpfs/ysm/project/fan/sb2723/00.database/hg38/Homo_sapiens.GRCh38.ncrna.fa
