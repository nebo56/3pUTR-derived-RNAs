'''
Created on Mar 16, 2014

@author: Nejc Haberman


Script will filter BED by lenght.
'''

import sys

def filterBEDbyCounts(fname_in, count_limit, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()

    while line:
        col = line.rstrip('\n').rsplit('\t')
        count = int(col[3])
        if abs(count) >= count_limit:
            fout.write(line)
        line = fin.readline()
    fout.close()
    fin.close()


if sys.argv.__len__() == 4:
    fname_in = sys.argv[1]
    count_limit = int(sys.argv[2])
    fname_out = sys.argv[3]
    filterBEDbyCounts(fname_in, count_limit, fname_out)
else:
    print("python filterBEDbyLenght.py <input_file> count_limit <output_file>")