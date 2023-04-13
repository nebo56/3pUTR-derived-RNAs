#!/bin/bash

target=$1 	# CAGE peaks in BED format
genome="./GRCh28.fa"
path="./"
motif="./motif-GGG.txt"

# flank 75 bps upstream and downstream from the peak CAGE position
python2 ${path}scripts/flankBEDpositionsCustom.py ${path}${target} ${path}${target}.flanked100.bed 100 100

# get fasta sequence
bedtools getfasta -s -fi ${genome} -fo ${path}${target}.flanked100.fa -bed ${path}${target}.flanked100.bed

# get coverage 
python ./scripts/miRNA-seed_coverage-norm-as_track.py ${path}${target}.flanked100.fa ${path}${motif} ${path}${target}.flanked100.GGG_enrichment.txt

# plot GGG enricment
Rscript ./scripts/plot_GGG_enruchment.R ${path}${target}.flanked100.GGG_enrichment.txt ${path}${target}.flanked100.GGG_enrichment.pdf