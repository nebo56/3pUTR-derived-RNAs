'''
Created on May 29, 2021

@author: Nejc Haberman

'''

import sys
import operator

# script will read the output from RNAfold and sum all '(' ')' as +1 normalise by the total number of positions
def sumRNAfold (fin_name):
    fin = open(fin_name)
    line = fin.readline()
    sum = None
    counter = 0
    while line:
        if line[0] == '>':
            counter += 1
            seq = fin.readline()
            fold = fin.readline()
            tokens = fold.rsplit(' ')
            fold = tokens[0]
            fold_str = fold.replace('.',"0,").replace('(',"1,").replace(')',"1,")   # replace '()' to 1, and '.' to 0
            fold_str = fold_str.rstrip(',')
            fold_str = fold_str.split(',')
            add = map(int, fold_str)
            if sum == None: #first time
                sum = add   #initialize the sum with the first array and convert it to int
            else:
                sum = map(operator.add, sum, add) #add new arry to sum
        line = fin.readline()
    fin.close()
    norm_sum = [x / float(counter) for x in sum]    #normalise by the number of positions
    print str(norm_sum).replace('[','').replace(']','').replace(' ','')

if sys.argv.__len__() == 2:
    fin = sys.argv[1]
    sumRNAfold(fin)
else:
    print "you need 1 arguments to run the script: input_fasta\n"
    quit()


