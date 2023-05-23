#!/bin/bash
#SBATCH -c 8
#SBATCH -J /mnt/storage/nejc/working/ENCODE-CAGE-K562/bowtie2-mapping.log
#SBATCH --mem 48G

# for nanopore sequencing analyses

bam1=$1 	# mapped BAM file of CAGE sample (replicate 1)
bam2=$2 	# mapped BAM file of CAGE sample (replicate 2)

genome="GRCh38"
thread="8"
path="./"

# convert BAM to BED
bedtools bamtobed -i ${path}${bam1} > ${path}${bam1}.bed
bedtools bamtobed -i ${path}${bam2} > ${path}${bam2}.bed

# 5' read positions
python2.7 ${path}scripts/getStart-BED.py ${path}${bam1}.bed ${path}${bam1}.5p.bed
python2.7 ${path}scripts/getStart-BED.py ${path}${bam2}.bed ${path}${bam2}.5p.bed

# sort 5' read positions
sort -k1,1 -k2,2n -k6,6 ${path}${bam1}.5p.bed > ${path}${bam1}.5p.sorted.bed
sort -k1,1 -k2,2n -k6,6 ${path}${bam2}.5p.bed > ${path}${bam2}.5p.sorted.bed

# get sum of 5' reads per position per replicate and for merged samples
python2.7 ${path}scripts/BEDsum.one_count.py ${path}${bam1}.5p.bed ${path}${bam1}.5p.sum.bed
python2.7 ${path}scripts/BEDsum.one_count.py ${path}${bam2}.5p.bed ${path}${bam2}.5p.sum.bed

# remove single sites
python2.7 ${path}scripts/filterBEDbyCounts.py ${path}${bam1}.5p.sum.bed 2 ${path}${bam1}.5p.sum-2min.bed
python2.7 ${path}scripts/filterBEDbyCounts.py ${path}${bam2}.5p.sum.bed 2 ${path}${bam2}.5p.sum-2min.bed

# merge replicates for clustering
cat ${path}${bam1}.5p.sum-2min.bed ${path}${bam2}.5p.sum-2min.bed | sort -k1,1 -k2,2n -k6,6 > ${path}${bam1}.${bam2}.merged.5p-2min.bed

# sum from both replicates after singletone removal
python2.7 ${path}scripts/BEDsum.py ${path}${bam1}.${bam2}.merged.5p-2min.bed ${path}${bam1}.${bam2}.merged.5p-2min.sum.bed

# get clusters in 20 nt window
bedtools merge -s -d 20 -c 6 -o distinct -i ${path}${bam1}.${bam2}.merged.5p-2min.sum.bed | awk '{print $1 "\t" $2 "\t" $3 "\t.\t.\t" $4}' > ${path}${bam1}.${bam2}.merged.5p-2min-merged20nt.bed

# find max CAGE 5' read peak per cluster
bedtools intersect -s -a ${path}${bam1}.${bam2}.merged.5p-2min.sum.bed -b ${path}${bam1}.${bam2}.merged.5p-2min-merged20nt.bed -wb > ${path}${bam1}.${bam2}.merged.5p-2min-merged20nt-peaks.bed
Rscript ${path}/scripts/get_max_cDNA_peak-inside_cluster.R ${path}${bam1}.${bam2}.merged.5p-2min-merged20nt-peaks.bed ${path}${bam1}.${bam2}.merged.5p-2min-merged20nt-max_peaks.bed
