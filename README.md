# DiBT-Notes
Learning notes for DBiT-Notes

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
        - Make shortcuts for the raw data folder if necessary, the folder name will be the same as hC2 (see examples)
        ```
        ln -s /gpfs/ysm/project/fan/sb2723/01.Spatial_hCortex/01.rawdata/usftp21.novogene.com/raw_data/hC2
        ```
        - Uzip and zip the raw data from Novogene using _gzip.sh_
   3. Set up the environment on HPC
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
   4. Make the index of the reference
   5. hh
   6. 


## Other stuffs
### HPC Work Commands
   1. Request for interactive job memory
   ```
   srun --pty -p interactive -c 1 --mem=6g bash
   ```
