#!/bin/bash

BED=$1 	# CAGE peaks in BED format (replicate 1)
genome="./GRCh38.fa"

# flank peak region to 100 bps upstream and downstream
python2 ./scripts/flankBEDpositionsCustom.py ${BED} ${BED}.flanked100.bed 100 100

# get fasta sequence
bedtools getfasta -s -fi ${genome} -bed ${BED}.flanked100.bed -fo ${BED}.flanked100.fa

# create results table
mkdir RESULTS.G4Hunter

# run G4Hunter prediction
python2 ./scripts/G4Hunter.py -i ${BED}.flanked100.fa -o RESULTS.G4Hunter -w 20 -s 1.2 

# in order to get the right format you need to remove specific lines
# $cat CAGE-Hs-K562-ENCSR000CJN-rep1-rep2-W25-S1.2.txt | grep -v "Start" | grep -v '>' > test.txt
# and add header "Start 	 End 	 Sequence	 Length 	 Score"

# plot 
Rscript ./scripts/plot_G4Hunter_score_as_metaplot.R ./RESULTS.G4Hunter/${BED}.flanked100-W20-S1.2.txt ./RESULTS.G4Hunter/${BED}.flanked100-W20-S1.2.pdf