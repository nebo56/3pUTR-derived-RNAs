#!/bin/bash -l
position=$1 #genomic positions in BED format
genome=./GRCh38.fa

# flank 75 bps in both directions
python ./scripts/flankBEDpositionsCustom.py ${position} ${position}.flanked75.bed 75 75

# get fasta
bedtools getfasta -s -fi ${genome} -bed ${position}.flanked75.bed -fo ${position}.flanked75.fasta

# get structure from RNAfold
RNAfold --noPS -i ${position}.flanked75.fasta > ${position}.flanked75.RNAfold.fasta

# sum of pairing probability
python ./scripts/sumRNAfold.py ${position}.flanked75.RNAfold.fasta > ${position}-RNAfold.sum.txt

# plot
Rscript ./scripts/plot_structure_enrichment.R ${position}-RNAfold.sum.txt ${position}-RNAfold.sum.pdf

#clean
rm ${position}.flanked75.fasta ${position}.flanked75.bed
