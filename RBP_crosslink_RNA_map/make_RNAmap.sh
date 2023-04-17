#!/bin/bash -l
# RNA-map of crosslink positions around the region. 

RBP=$1	# RBP crosslink positions in BED format (each crosslink needs to be in a newe line)
target_positions=$2	# target positions in BED format
name=$3

mkdir ${name}

# flank the region of your interest for 100 bps upstream and downstream
python ./scripts/flankBEDpositionsCustom.py ${target_positions} ${target_positions}.flanked100.bed 100 100

# get crosslink coverage coverage (if sorted use '-sorted' option)
bedtools coverage -s -b ${clusters} -a ${target_positions}.flanked100.bed -d > ${target_positions}.flanked100.CL-coverage.bed

# plot RNA map 
Rscript ./scripts/draw-RNA-maps.R ${target_positions}.flanked100.CL-coverage.bed ${name}

# move results into new folder
mv ${name}.pdf ${name}

# clean
rm ${target_positions}.flanked100.bed
rm ${target_positions}.flanked100.CL-coverage.bed
