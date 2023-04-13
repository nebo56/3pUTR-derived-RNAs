#!/bin/bash -l

# set correct paths for the following tools
export PATH=~/Programs/bedtools2.22.1/bin:$PATH
export PATH=~/Programs/samtools-0.1.19:$PATH
export PATH=~/Programs/scripts:$PATH
export PATH=~/Programs/fastx_toolkit_0.0.13:$PATH

data1=$1        # compressed FASTQ - first pair
data2=$2        # compressed FASTQ - second pair
path=`pwd -P`
path=${path}/
genome=./GRCh38.fa
index=./GRCh38.fa.fai
chr_GTF="./gencode.v27.primary_assembly.annotation.gtf"
thread="8"

mkdir $path$data1$data2-STAR

# adapter removal
cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o $path$data1-adapter_trimmed.fastq.gz -p $path$data2-adapter_trimmed.fastq.gz $path$data1 $path$data2 > $path$data1$data2.adapterTrim.metrics

# second adapter removal - double ligation 
cutadapt --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o $path$data1-adapter_trimmed.round2.fastq.gz -p $path$data2-adapter_trimmed.round2.fastq.gz $path$data1-adapter_trimmed.fastq.gz $path$data2-adapter_trimmed.fastq.gz > $path$data1$data2.adapterTrim.round2.metrics

# mapping to genome
time STAR --runMode alignReads --runThreadN $thread --genomeDir $genome_dir --readFilesCommand zcat --readFilesIn $path$data1-adapter_trimmed.round2.fastq.gz $path$data2-adapter_trimmed.round2.fastq.gz --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix $path$data1$data2-STAR/ --outSAMattributes All --outStd BAM_SortedByCoordinate --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo

# SAM to BAM
samtools view -hSb $path$data1$data2-STAR/Aligned.out.sam > $path$data1$data2-STAR/$data1.$data2.bam
rm $path$data1$data2-STAR/Aligned.out.sam

# PCR duplicate removal
python2 ${path}scripts/barcode_collapse_pe.py --bam $path$data1$data2-STAR/$data1.$data2.bam --out_file $path$data1$data2-STAR/$data1.$data2.rmDup.bam --metrics_file $path$data1.$data2-STAR/$data1.$data2.rmDup.metrics

# remove low quality reads
samtools view -hb -q255 $path$data1$data2-STAR/$data1.$data2.rmDup.bam > $path$data1$data2-STAR/$data1.$data2.rmDup-q255.bam

# BAM to BED
bedtools bamtobed -bedpe -mate1 -i $path$data1$data2-STAR/$data1.$data2.rmDup-q255.bam | sort -k1,1 -k2,2n -k6,6 > $path$data1$data2-STAR/$data1.$data2.rmDup-q255.bed

# BAM to single-end BED format
python ${path}scripts/pairedBED2singleBED.py $path$data1$data2-STAR/$data1.$data2.rmDup-q255.bed $path$data1$data2-STAR/$data1.$data2.rmDup-q255-single_end.bed

# get cross link positions and convert it to bedGraph
python ./scripts/BEDtoXlink.py $path$data1$data2-STAR/$data1.$data2.rmDup-q255-single_end.bed ${path}${data}-xlink.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-xlink.bed > ${path}${data}-xlink-sorted.bed
python ./scripts/BEDsum.py ${path}${data}-xlink-sorted.bed ${path}${data}-xlink-sorted-sum.bed
python ./scripts/BED2BEDgraph.py ${path}${data}-xlink-sorted-sum.bed ${path}${data}-xlink-sorted-sum.bedGraph
rm ${path}${data}-xlink.bed

# print unique number of cDNAs
echo "Unique number of cDNAs" > ${path}${data}-01.Mapping-REPORT.txt
wc -l ${path}${data} >> ${path}${data}-01.Mapping-REPORT.txt

# compress the original data
gzip ${path}${data}
gzip ${path}${data}-xlink-sorted-sum.bed
gzip ${path}${data}-xlink-sorted-sum.bedGraph
