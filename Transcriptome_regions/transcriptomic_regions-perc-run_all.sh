#!/bin/bash

# script will return percentage of reads per genomic region (3'UTR, 5'UTR, CDS, intron)
# INPUT folder = FILE_PATH: sorted crosslink BED file, where each line is one read

OUTPUT=$1	#text file for results
FILE_PATH=/Users/nhaberma/Imperial/CAGE-genomic_regions-all_cell_lines-unique_positions-pipeline/*-uniq.bed
gen_regions=gencode.v26.annotation.seg.CDS.UTR3.UTR5.intron.merged.bed

echo  "sample	read_number	perc.3UTR	perc.5UTR	perc.CDS	perc.intron" >> $OUTPUT

for single_position in $FILE_PATH
do
  bedtools intersect -s -a $single_position -b $gen_regions -wb | uniq > ${single_position}-regions.bed
  nsingle_position=`wc -l < "${single_position}-regions.bed"`
  utr3=`grep -c UTR3 "${single_position}-regions.bed"`
  utr5=`grep -c UTR5 "${single_position}-regions.bed"`
  cds=`grep -c CDS "${single_position}-regions.bed"`
  intron=`grep -c intron "${single_position}-regions.bed"`
  utr3perc=`echo "scale=2; $utr3*100/$nsingle_position" | bc`
  utr5perc=`echo "scale=2; $utr5*100/$nsingle_position" | bc`
  cds_perc=`echo "scale=2; $cds*100/$nsingle_position" | bc`
  intron_perc=`echo "scale=2; $intron*100/$nsingle_position" | bc`

  filename=${single_position##*/}     #remove the path from the filename
  RBP=$(echo $filename| cut -d'.' -f 1) #get sample's name out of the file name vy removing everything after "_"
  
  echo $RBP'	'$nsingle_position'	'$utr3perc'	'$utr5perc'	'$cds_perc'	'$intron_perc >> $OUTPUT
  gzip $single_position
  rm ${single_position}-regions.bed
done

