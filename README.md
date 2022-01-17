# DiBT-Notes
Learning notes for DBiT-Notes

## DBiT work flow
  1. Experiments on tissue samples;
     - Sample thickness < 10 um;
     - 50 * 50 barcodes;
     - Extract the Sequences using Illumina kit;
    
   2. Sequencing by Novogene;
     - Check the library QC report.
     - Batch download the data to HPC folder;
       - Make new project folder and put raw data in it with a name of 0.raw_data
      ```
      wget -r -c ftpxxx
      ```
   3. Set up the environment on HPC
      - Follow the instruction of ST Pipeline https://github.com/jfnavarro/st_pipeline
