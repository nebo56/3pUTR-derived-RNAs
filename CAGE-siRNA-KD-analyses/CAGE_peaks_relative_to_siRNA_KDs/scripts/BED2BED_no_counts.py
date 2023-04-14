'''
Created on Nov 13, 2013

@author: Nejc Haberman

'''


import sys

def BED2BED_no_counts(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        count = int(col[4])
        
        for i in range(0,int(count)):
            fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + "" + '\t' + col[4] + "\t" + col[5] + '\n')
        line = fin.readline()
                    
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BED2BED_no_counts(fname_in, fname_out)
else:
    print("python BED2BED_no_counts.py <input_file> <output_file>")