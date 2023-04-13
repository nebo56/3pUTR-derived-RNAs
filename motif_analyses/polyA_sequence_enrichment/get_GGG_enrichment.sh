#!/bin/bash

target=$1 	# CAGE peaks in BED format
genome="./GRCh28.fa"
path="./"
motif="./motif-AAUAAA.tab"

# get fasta sequence
bedtools getfasta -s -fi ${genome} -fo CAGE-Hs-K562-ENCSR000CJN-rep1-rep2.max_peak.filtered.3UTR-flanked100.fa -bed ./data/CAGE-Hs-K562-ENCSR000CJN-rep1-rep2.max_peak.filtered.3UTR-flanked100.bed
bedtools getfasta -s -fi ${genome} -fo gencode.v27.basic.3UTR-End-flanked100.fa -bed ./data/gencode.v27.basic.3UTR-End-flanked100.bed

# get coverage 
python2 ./scripts/k-mer_coverage-norm.py CAGE-Hs-K562-ENCSR000CJN-rep1-rep2.max_peak.filtered.3UTR-flanked100.fa ${path}${motif} CAGE-Hs-K562-ENCSR000CJN-rep1-rep2.max_peak.filtered.3UTR-flanked100.polyA.enrichment.txt
python2 ./scripts/k-mer_coverage-norm.py gencode.v27.basic.3UTR-End-flanked100.fa ${path}${motif} gencode.v27.basic.3UTR-End-flanked100.polyA.enrichment.txt

# plot polyA motif enricment
Rscript ./scripts/motif-enrichment-polyA-3UTR-CAGE-K562_max_peak_3pUTR-end.R ${path}${target}.flanked100.GGG_enrichment.txt ${path}${target}.flanked100.GGG_enrichment.pdf