sample=$1
indir=$2
FW=$indir/$sample/$sample.R1.fq.gz
RV=$indir/$sample/$sample.R2.fq.gz
MAP=/gpfs/ysm/project/fan/sb2723/00.database/hg38/STARindex
ANN=/gpfs/ysm/project/fan/sb2723/00.database/hg38/gencode.v39.annotation.gtf
CONT=/gpfs/ysm/project/fan/sb2723/00.database/hg38/STARindex_nc
ID=/gpfs/ysm/project/fan/sb2723/00.database/barcodes-AB.xls
OUTPUT=$indir/03.stpipeline/$sample
mkdir -p $indir/03.stpipeline/$sample
TMP=$indir/03.stpipeline/$sample/tmp
mkdir -p $indir/03.stpipeline/$sample/tmp
EXP=$sample

st_pipeline_run.py \
  --output-folder $OUTPUT \
  --temp-folder $TMP \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --ids $ID \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $EXP \
  --htseq-no-ambiguous \
  --verbose \
  --threads 16 \
  --log-file $OUTPUT/${EXP}_log.txt \
  --star-two-pass-mode \
  --no-clean-up \
  --contaminant-index $CONT \
  --disable-clipping \
  --min-length-qual-trimming 30 \
  $FW $RV
