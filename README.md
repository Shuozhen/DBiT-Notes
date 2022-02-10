# DiBT-Notes
Learning notes for DBiT-Notes, credit to Mingyu Yang https://github.com/MingyuYang-Yale/DBiT-seq.

## DBiT work flow
  1. Experiments on tissue samples;
     - Sample thickness < 10 um;
     - 50 * 50 barcodes;
     - Extract the Sequences using Illumina kit;
  2. Sequencing by Novogene;
     - Check the library QC report.
       - Copy the ![image](https://user-images.githubusercontent.com/25277637/152657283-475b4575-fe49-41ac-8979-a2cfa375c31b.png) into the excel file named "Sample QC Overall"
     - Batch download the data to HPC folder;
       - Make new project folder and put raw data in it with a name of **0.raw_data**
        ```
        wget -r -c ftpxxx
        ```
       - Make shortcuts for the raw data folder if necessary, the folder name will be the same as hC2
        ```
        ln -s /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/01.rawdata/usftp21.novogene.com/raw_data/hC2
        ```
       - Uzip and zip the raw data from Novogene using _gzip.sh_
   3. Set up the ST Pipeline environment (conda) on HPC
      - Install X11
        ```
        launchctl load -w /Library/LaunchAgents/org.macosforge.xquartz.startx.plist
        echo $DISPLAY
        ```
      - Follow the instruction of ST Pipeline https://github.com/jfnavarro/st_pipeline
      - Current procedure
        ```
        module load miniconda
        conda create -n st-pipeline python=3.7
        conda activate st-pipeline
        conda install PySam
        conda install Numpy
        conda install Cython
        pip install taggd
        pip install stpipeline
        ```
      - Test the installment
        ```
        st_pipeline_run.py -h
        ```
   4. Set up Perl environment on HPC, and run _effective.sh_ afterwards
        ```
        wget https://cpan.metacpan.org/authors/id/N/NW/NWCLARK/PerlIO-gzip-0.20.tar.gz
        module avail Perl
        module load Perl/5.28.0-GCCcore-7.3.0
        tar -zxvf PerlIO-gzip-0.20.tar.gz 
        cd PerlIO-gzip-0.20/
        mkdir mybuild
        perl Makefile.PL PREFIX=/gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/00.bin/PerlIO-gzip-0.20/mybuild
        make
        make install
        ```
## HPC Data Processing
   1. Make the index of the reference
      - Once it's settled up, no need to change unless there's new sample or new updates
        - Current version of human gene reference
        ```
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
        wget http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.long_noncoding_RNAs.gtf.gz
        wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
        ```
        - install the STAR and samtool to establish the reference database, better in a database folder
        ```
        module load miniconda
        conda activate st-pipeline
        conda install -c bioconda star
        conda install -c bioconda samtools openssl=1.0
        st_pipeline_run.py -v
        rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ .
        ```
        - Unzip the annotation file and combine them into one, delete the original seperate files (not double confirmed)
        ```
        for i in {1..22} X Y M; do gzip -d chr$i.fa.gz;done
        for i in {1..22} X Y M; do cat chr$i.fa; done >> hg38.fa
        for i in {1..22} X Y M; do rm chr$i.fa;done
        cut -f1  Homo_sapiens.GRCh38.105.gtf | uniq
        cut -f1 gencode.v39.annotation.gtf | uniq
        grep '>chr' hg38.fa
        ```
        - Make the directory of STARindex and STARindex_nc under the genome directory
        - The following in the _starindex.sh_ and _starindex_nc.sh_ needed to be changed, and the sbatch both files, remember to request for the interactive nodes for the job
          - Make a directory of STARindex and set it as the genomeDir
          ```
          --genomeDir /gpfs/ysm/project/fan/sb2723/00.database/hg38/STARindex
          ```
          - The genome fasta file is the one we combined last step
          ```
           --genomeFastaFiles /gpfs/ysm/project/fan/sb2723/00.database/hg38/hg38.fa
          ```
          - The annotation file should be coordinate with the genome fasta file
            - The chromosome name should be chr* or numbers/X/Y only
          ```
          --sjdbGTFfile /gpfs/ysm/project/fan/sb2723/00.database/hg38/gencode.v39.annotation.gtf 
          ```
          - This line denotes the sequencing length, we're doing 150 here
          ```
          --sjdbOverhang 149
          ```
          - The limit of genome generate RAM should be adjusted by the instruction (* double confirm with Mingyu)
          ```
          --limitGenomeGenerateRAM 50000000000
          ```
          - The directory should be changed to nc folder and the fasta files should also be changed.
          ```
          --genomeDir /gpfs/ysm/project/fan/sb2723/00.database/hg38/StarIndex_nc
          --genomeFastaFiles /gpfs/ysm/project/fan/sb2723/00.database/hg38/Homo_sapiens.GRCh38.ncrna.fa
          ```
   2. Filter the raw data and rearrange read format to be compatible with ST Pipeline using _effective.sh_
      - Perl file is used for the processing, _1-effective.pl_
        ```
        sbatch hC1.stpipeline.sh
        ```
      - Cat the log file to paste the barcode A and B and barcode AB into the excel file named "Sample QC Overall"
        - Calculate the Effective reads by barcodeAB/(Raw Reads/s)*100%
   3. Run st-pipeline
      - Remember to get into the miniconda environment
        ```
        module load miniconda
        conda activate st-pipeline
        ```
      - Run ST pipeline.
        - _stpipeline.sh_ is in /00.sh, just run _hC1.stpipeline.sh_ because it's running the _stpipeline.sh_
        ```
        sbatch hC1.stpipeline.sh
        ```
   4. Change ID using _changeid.sh_
      - Run the _changeid.sh_ for samples.
        ```
        module load miniconda
        conda activate st-pipeline
        srun --pty -p interactive --mem=16g bash
        sh changeid.sh
        ```
## Raw Image Processing 
https://github.com/edicliuyang/DBiT-seq_FFPE/blob/master/Figure_Processing/Pixel_identification.m
   1. Crop the image using PS
   2. Use Image -> Adjustment -> Threshold
   3. Use Image -> Adjustment -> Invert
   4. Use the Matlab script _Pixel_identification.m_ to generate position information

## R scripts


## Other stuffs
### HPC Work Commands
   1. Request for interactive job memory
   
       ```
       srun --pty -p interactive -c 1 --mem=6g bash
       ```
   2. search in the history
     
       ```
       history | grep interactive
       ```
   3. Substitute something in batch
     
       ```
       :%s/(original)/(replaced)/g
       ```
   4. Save history
     
       ```
       history > logxxx
       ```
   5. Change the file permission
      ```
      chmod 755 -R STARindex_nc
      ```
   6. Environment setting and shortcuts
      ```
      vi ./~bashrc
      ```
