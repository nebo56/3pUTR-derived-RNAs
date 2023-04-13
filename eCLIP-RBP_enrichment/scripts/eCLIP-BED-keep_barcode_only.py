'''
Created on Oct 28, 2019

@author: Nejc Haberman

Script will remove info header from the 4th column and kept only the barcode
example: bedtools bamtobed -i eCLIP-K562-FAM120A-rep1_1.fastq.bam > eCLIP-K562-FAM120A-rep1_1.fastq.bed

- before:
chr1    88837882        88837916        TTCCG:HWI-D00611:119:C6K7PANXX:2:1314:18227:81523       255     -

- after:
chr1    88837882        88837916        TTCCG       255     -
'''


import sys

def filterBED(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    #header = fin.readline()
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        info = col[3].rsplit(':')
        barcode = info[0]
        fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + barcode + '\t' + col[4] + '\t' + col[5] + '\n')
        line = fin.readline()
                    
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    filterBED(fname_in, fname_out)
else:
    print("python filterBED.py <input_file> <output_file>")
