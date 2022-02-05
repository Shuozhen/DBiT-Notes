#!/bin/bash
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuozhen.bao@yale.edu
#SBATCH --job-name=test
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=4g 
#SBATCH --time=120:00:00

module load Perl/5.28.0-GCCcore-7.3.0
export PERL5LIB=/gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/00.bin/PerlIO-gzip-0.20/mybuild/lib/perl5/site_perl/5.28.0/x86_64-linux-thread-multi:$PERL5LIB

#perl /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/00.bin/1-effective.pl  -indir /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/01.rawdata -outdir /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/02.effective -insertsize hC1

perl /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/00.bin/1-effective.pl  -indir /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/01.rawdata -outdir /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/02.effective -insertsize hC2
