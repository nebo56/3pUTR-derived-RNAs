#!/bin/bash
# weblogo tools needs to be installed first from https://weblogo.berkeley.edu/

position=$1 #genomic positions of CAGE peaks in BED format
genome=./GRCh38.fa

# flank region for 5 bps upstream and downstream
python2 ./scripts/flankBEDpositionsCustom.py ${1} ${1}-flanked5.bed 5 5

# get fasta sequence around CAGE peaks
bedtools getfasta -s -fi ~/UCL/genome-data/GRCh38.p10.genome.fa -bed ${1}-flanked5.bed -fo ${1}-flanked5.fasta

# get web logos
~/Programs/weblogo/seqlogo -S -Y -c -n -s -5 -f ${1}-flanked5.fasta > ${1}-flanked5.eps