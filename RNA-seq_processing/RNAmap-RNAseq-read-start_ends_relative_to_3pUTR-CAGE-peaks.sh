#!/bin/bash

# script will return coverage of 5' and 3' end RNA-seq positions around 3'UTR CAGE peaks
# run it for each replicate seperately

cage=$1 	# 3'UTR CAGE peaks in BED format
rnaseq_readStarts=$2 	# 5'ends of RNAseq reads
rnaseq_readEnd=$3 	# 3'ends of RNAseq reads
path="./"

# expand the targeting region to 75 bps upstream and downstream
python2 ${path}scripts/flankBEDpositionsCustom.py ${path}${cage} 75 75 ${path}${cage}.flanked75.bed

# get coverage of RNA-seq read-starts and read-ends around CAGE peaks
bedtools coverage -sorted -s -b {path}${rnaseq_readEnds} -a ${path}${cage}.flanked75.bed -d > ${path}${cage}.flanked75.RNAseq-read_ends.bed
bedtools coverage -sorted -s -b {path}${rnaseq_readStarts} -a ${path}${cage}.flanked75.bed -d > ${path}${cage}.flanked75.RNAseq-read_starts.bed

# plot RNAseq read-start and read-end positions relative to 3'UTR CAGE peaks
Rscript ${path}scripts/Plot_RNAseq_reads_relative_to_CAGE_peaks.R ${path}${cage}.flanked75.RNAseq-read_ends.bed ${path}${cage}.flanked75.RNAseq-read_starts.bed RNAseq-read_starts_ends_relative_to_CAGE_peaks.pdf
