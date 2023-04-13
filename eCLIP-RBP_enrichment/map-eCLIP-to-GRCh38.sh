#!/bin/bash -l

# set correct paths for the following tools
export PATH=~/Programs/bedtools2.22.1/bin:$PATH
export PATH=~/Programs/samtools-0.1.19:$PATH
export PATH=~/Programs/scripts:$PATH
export PATH=~/Programs/fastx_toolkit_0.0.13:$PATH

data=$1
path=`pwd -P`
path=${path}/
genome=./GRCh38.fa

gunzip ${path}${data}.gz

# get cross link positions and convert it to bedGraph
python ./scripts/BEDtoXlink-CLIP.py ${path}${data} ${path}${data}-xlink.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-xlink.bed > ${path}${data}-xlink-sorted.bed
python ./scripts/BEDsum.py ${path}${data}-xlink-sorted.bed ${path}${data}-xlink-sorted-sum.bed
python ./scripts/BED2BEDgraph.py ${path}${data}-xlink-sorted-sum.bed ${path}${data}-xlink-sorted-sum.bedGraph
rm ${path}${data}-xlink.bed

# print unique number of cDNAs
echo "Unique number of cDNAs" > ${path}${data}-01.Mapping-REPORT.txt
wc -l ${path}${data} >> ${path}${data}-01.Mapping-REPORT.txt

# count cDNAs that have more then 1 crosslink (cDNA start > 1)
python ./scripts/count_cDNAs.py ${path}${data}-xlink-sorted-sum.bed >> ${path}${data}-01.Mapping-REPORT.txt

# compress the original data
gzip ${path}${data}
gzip ${path}${data}-xlink-sorted-sum.bed
gzip ${path}${data}-xlink-sorted-sum.bedGraph
