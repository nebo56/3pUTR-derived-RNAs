'''
Created on Nov 19, 2019
 
@author: Nejc Haberman
 
description:
Script will accept fasta sequnces (all sequences must have the same length or at least the same begining) 
where we look for the coverage based on all motifs (each motif must be in a new line).
For each sequnce we will look for the coverage (each nt overlapping with motif will count as 1) 
and sum together coverage from all fasta sequences and divided by number of sequences and their length as an
average distribution of kmers

!!!! each kmer will be done separately and exported into separate file !!!
 
!!!! here motif file contains 2 columns: motif, miRNA name

!!!! we clip 7mer to 6 mer from both sides and select the one that gives us a better coverage

input:
- input_fasta
- motifs (kmers)
- output_fasta
 
output:
- fasta with extra line of coverage (nt) per sequence
- total coverage 1 nt resolution of all fasta files 
'''

import sys
 
# load motifs from the file !!!! here motif file contains 3 columns: motif, strand, miRNA name
def load_motifs(fin_fname_motifs):
    fin = open(fin_fname_motifs, "rt")
    motifs = []
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        motif = str(col[0]).upper()
        name = col[1]
        new = (motif, name)
        motifs.append(new)
        line = fin.readline()
    fin.close()
    return motifs
 
# set the coverage of the sequence
def set_coverage(seq, motif, coverage):
    length = motif.__len__()
    motif_pos = seq.find(motif)
    pos = motif_pos
    while motif_pos != -1:
        coverage[pos:pos+length] = length * [1]
        motif_pos = seq[pos+1:].find(motif) #we search for the same motif downstream fro mthe previous one
        pos = motif_pos + pos + 1
    return coverage

def sum_coverage(total_coverage, adding_coverage):
    for i in range(0, total_coverage.__len__()):
        total_coverage[i] = total_coverage[i] + adding_coverage[i]
    return total_coverage 
 
def get_coverage(fin_fname, fin_fname_motifs, fout_fname):
    motifs = load_motifs(fin_fname_motifs)
    info = ""    
    for i in range(0,motifs.__len__()): #for each motif we do a new search and exported to a new file
        fin = open(fin_fname, "rt")
        new = motifs[i]
        motif = new[0]

        motif_1 = motif[1::]    # we trim the first nucleotide
        motif_2 = motif[:-1]    # we trim the last nucleotide

        name = new[1]
        fout = open(fout_fname + "." + motif + ".coverage.txt", "w")
        counter = 0
        coverage_sum_1 = 0
        coverage_sum_2 = 0
        total_coverage_1 = 0
        total_coverage_2 = 0
        seq_length = 0
        coverage_1 = [0] * seq_length     #initialize
        coverage_2 = [0] * seq_length  # initialize
        line = fin.readline()
        while line:
            if line[0] == '>':
                info = line
            else:
                counter += 1
                seq = line.rstrip('\n')
                seq = str(seq).upper()
                if seq_length == 0:    #first time we need to get the length of the sequence
                    seq_length = seq.__len__()
                    coverage_1 = [0] * seq_length
                    coverage_2 = [0] * seq_length
                    total_coverage_1 = [0] * seq_length
                    total_coverage_2 = [0] * seq_length
                set_coverage(seq, motif_1, coverage_1)
                set_coverage(seq, motif_2, coverage_2)
                total_coverage_1 = sum_coverage(total_coverage_1, coverage_1)
                total_coverage_2 = sum_coverage(total_coverage_2, coverage_2)
                coverage_1 = [0] * seq_length
                coverage_2 = [0] * seq_length
            line = fin.readline()

        if max(total_coverage_1) >= max(total_coverage_2):  # we select the trimmed seed that has a stronger peak of enrichment
            total_coverage = total_coverage_1
        else:
            total_coverage = total_coverage_2

        norm_total_coverage = [x / float(counter) for x in total_coverage]  #normalize by the number of sequences
        seed = str(new).replace("'","").replace("(","").replace(")","").replace(", ","\t")  #motif \t strand \t name
        fout.write(str(norm_total_coverage).replace(' ','').replace('[','').replace(']','').replace(",",str("\t" + seed + '\n')) + '\t' + seed + '\n')
        fin.close()
        fout.close()


'''
fin_fname_fasta = "/Users/nhaberma/Imperial/eiCLIP-Chris/Ago2-HeLa-miRNA-seed-coverage_mir-TargetScan-2/ams-a2-hela_trimmed_single-sum-TARGET.bed.flank30.fasta"
fin_fname_motifs = "/Users/nhaberma/Imperial/eiCLIP-Chris/Ago2-HeLa-miRNA-seed-coverage_mir-TargetScan-2/TargetScanHuman.Release.7.2-miR_Family_Info-HeLa-mean10reads.seed.rev_comp.name.tab"
fout_fname = "/Users/nhaberma/Imperial/eiCLIP-Chris/Ago2-HeLa-miRNA-seed-coverage_mir-TargetScan-2/test.txt"
get_coverage(fin_fname_fasta, fin_fname_motifs, fout_fname)
'''

if sys.argv.__len__() == 4:
    fin_fname_fasta = sys.argv[1]
    fin_fname_motifs = sys.argv[2]
    fout_fname = sys.argv[3]
    get_coverage(fin_fname_fasta, fin_fname_motifs, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t3 arguments are needed\n" + '\n' +"example:\t $ python k-mer_coverage.py input_fname.fasta motifs.tab output_fname.fasta"    


