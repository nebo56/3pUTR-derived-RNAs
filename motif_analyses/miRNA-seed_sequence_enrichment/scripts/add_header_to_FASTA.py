'''
Created on Aug 7, 2019

@author: Nejc

Script will add extra line as header ">" to each sequence.
'''

import sys

def add_header(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    line = fin.readline()
    while line:
        fout.write(">" + line)
        fout.write(line)
        line = fin.readline()
    fin.close()
    fout.close()

if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    add_header(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n"
