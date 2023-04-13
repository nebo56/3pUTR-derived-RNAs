'''
Created on Aug 20, 2013

@author: Nejc

The script will write end of the cluster and start positions to a different file
'''

import sys

def set_CrossLinks(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    foutMiddle = open(fout_fname + "-Middle.bed", "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        pos1 = col[1]
        pos2 = col[2]
        #cDNA = col[3]
        cDNA = ""
        strand = col[5]
        
        if strand == "+":
            start = int(pos1)
            end = int(pos2)
        elif strand == "-":
            start = int(pos2)
            end = int(pos1)
            
        middle = int(pos1) + int((int(pos2) - int(pos1)) / 2)
        foutMiddle.write(chr + '\t' + str(middle) + '\t' + str(middle+1) + '\t' + cDNA + '\t' + "" + '\t' + strand + '\n')
        line = fin.readline()
    foutMiddle.close()

if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    set_CrossLinks(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python get_cDNA_middle.py input_fname.bed output_fname"
