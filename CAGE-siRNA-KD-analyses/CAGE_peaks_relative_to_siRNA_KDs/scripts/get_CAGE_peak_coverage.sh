#!/bin/bash -l

target=$1	#si RNA start
CAGE1=$2		# ctss rep1
CAGE2=$3		# ctss rep2
CAGE3=$4		# ctss rep3
name=$5

# flank target sequence
python2 ./scripts/flankBEDpositionsCustom.py ${target} ${target}-flank100.bed 100 100

# get 5' read counts for replicate 1
bedtools intersect -s -a ${CAGE1} -b ${target}-flank100.bed -wa | uniq > ${CAGE1}-targets.bed
python2 ./scripts/BED2BED_no_counts.py ${CAGE1}-targets.bed ${CAGE1}-targets-no_counts.bed

# get 5' read counts for replicate 2
bedtools intersect -s -a ${CAGE2} -b ${target}-flank100.bed -wa | uniq > ${CAGE2}-targets.bed
python2 ./scripts/BED2BED_no_counts.py ${CAGE2}-targets.bed ${CAGE2}-targets-no_counts.bed

# get 5' read counts for replicate 3
bedtools intersect -s -a ${CAGE3} -b ${target}-flank100.bed -wa | uniq > ${CAGE3}-targets.bed
python2 ./scripts/BED2BED_no_counts.py ${CAGE3}-targets.bed ${CAGE3}-targets-no_counts.bed

# get coverage around the target region
bedtools coverage -s -a ${target}-flank100.bed -b ${CAGE1}-targets-no_counts.bed -d > ${target}-flank100-${CAGE1}.coverage.bed
bedtools coverage -s -a ${target}-flank100.bed -b ${CAGE2}-targets-no_counts.bed -d > ${target}-flank100-${CAGE2}.coverage.bed
bedtools coverage -s -a ${target}-flank100.bed -b ${CAGE3}-targets-no_counts.bed -d > ${target}-flank100-${CAGE3}.coverage.bed

# control
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1A_rep1.ctss.bed -d > ${target}-flank100-control_batch1A_rep1.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1A_rep2.ctss.bed -d > ${target}-flank100-control_batch1A_rep2.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1A_rep3.ctss.bed -d > ${target}-flank100-control_batch1A_rep3.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1B_rep1.ctss.bed -d > ${target}-flank100-control_batch1B_rep1.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1B_rep2.ctss.bed -d > ${target}-flank100-control_batch1B_rep2.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch1B_rep3.ctss.bed -d > ${target}-flank100-control_batch1B_rep3.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch2_rep1.ctss.bed -d > ${target}-flank100-control_batch2_rep1.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch2_rep2.ctss.bed -d > ${target}-flank100-control_batch2_rep2.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch2_rep3.ctss.bed -d > ${target}-flank100-control_batch2_rep3.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch3_rep1.ctss.bed -d > ${target}-flank100-control_batch3_rep1.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch3_rep2.ctss.bed -d > ${target}-flank100-control_batch3_rep2.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch3_rep3.ctss.bed -d > ${target}-flank100-control_batch3_rep3.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch4_rep1.ctss.bed -d > ${target}-flank100-control_batch4_rep1.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch4_rep2.ctss.bed -d > ${target}-flank100-control_batch4_rep2.bed
bedtools coverage -s -a ${target}-flank100.bed -b Negative_control_batch4_rep3.ctss.bed -d > ${target}-flank100-control_batch4_rep3.bed

Rscript ./scripts/draw_RNAmap.R ${target}-flank100-${CAGE1}.coverage.bed ${target}-flank100-${CAGE2}.coverage.bed ${target}-flank100-${CAGE3}.coverage.bed ${name}

rm ${target}-flank100.bed 
rm ${CAGE1}-targets.bed ${CAGE1}-targets-no_counts.bed
rm ${CAGE2}-targets.bed ${CAGE2}-targets-no_counts.bed
rm ${CAGE3}-targets.bed ${CAGE3}-targets-no_counts.bed
rm ${CAGE4}-targets.bed ${CAGE4}-targets-no_counts.bed
