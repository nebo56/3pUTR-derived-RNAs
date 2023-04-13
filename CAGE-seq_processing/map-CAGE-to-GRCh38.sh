#!/bin/bash

# genome annotation available here https://www.gencodegenes.org/human/

fastq=$1 	#compressed FASTQ
path=./
genome=./GRCh38.p10.genome.fa
chr_GTF="./gencode.v27.primary_assembly.annotation.gtf"
genome="GRCh38"
thread="8"
genome_dir="./GRCh38.p10-STAR"

mkdir $path$fastq-STAR

# map CAGE FASTq reads to GRCh38 genome using STAR alignment and not allowing trimming on 5' end of the read
time STAR --runMode alignReads --runThreadN $thread --readFilesCommand zcat --genomeDir $genome_dir --readFilesIn ${path}${fastq} --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix $path$fastq-STAR/ --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType Extend5pOfRead1 --clip5pNbases 9

# move mapped BAM file into the folder
mv $path$fastq-STAR/Aligned.sortedByCoord.out.bam $path$fastq-STAR/$fastq.bam

# create index
samtools index $path$fastq-STAR/$fastq.bam $path$fastq-STAR/$fastq.bam.bai

