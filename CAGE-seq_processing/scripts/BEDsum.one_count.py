'''
Created on Mar 16, 2014

@author: Nejc Haberman


Script will sum together BED but BED needs to be sorted and each line will be counted as 1.
'''

import sys

def BEDsum(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    lastChr = None
    lastStart = None
    lastEnd = None
    lastStrand = None
    cDNAsum = 0
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        cDNA = 1
        strand = col[5]
        if lastChr == None or (lastChr == chr and lastStart == start and lastEnd == end and lastStrand == strand):
            cDNAsum += int(cDNA)
        else:
            fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + str(cDNAsum) + '\t.\t' + lastStrand + '\n')
            cDNAsum = int(cDNA)
        lastChr = chr
        lastStart = start
        lastEnd = end
        lastStrand = strand
        line = fin.readline()
    fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + str(cDNAsum) + '\t.\t' + lastStrand + '\n')
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDsum(fname_in, fname_out)
else:
    print("python BEDsum.py <input_file> <output_file>")
