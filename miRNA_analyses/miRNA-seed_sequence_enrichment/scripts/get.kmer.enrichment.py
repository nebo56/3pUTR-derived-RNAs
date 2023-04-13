'''
Created on May 10, 2017

@author: Nejc Haberman

Script will calculate the kmer enrichment between each RBP and control kmer coverage. Control needs to be in the last column.

''' 

import sys

def get_enrichment (fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    header = fin.readline()
    new_header = header.replace("coverage","enrichment")
    fout.write(new_header)
    line = fin.readline()
    while line:
        col = line.rstrip("\n").rsplit("\t")
        fout.write(col[0])  #tetramer
        for i in range(1, col.__len__()-1):
            fout.write("\t" + str(float(col[i]) / float(col[-1])))
        fout.write("\n")
        line = fin.readline()
    fout.close()
    fin.close()

if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    get_enrichment(fname_in, fname_out)
else:
    print("python get.kmer.enrichment.py <input_file> <output_file>")

'''
add_RPKM("/media/skgthab/storage/UCL/2017.04.28-CLIPo-iONMF-ENCODE/CLIPo-eCLIP-K562-merged-narrow_clusters-RBPs_and_mock/RBPs/RPKMs/eCLIP-K562-all_RBPs_and_mocks-max_peaks-narrow-merged-coverage_intersect-narrow_clusters.tab","/media/skgthab/storage/UCL/2017.04.28-CLIPo-iONMF-ENCODE/CLIPo-eCLIP-K562-merged-narrow_clusters-RBPs_and_mock/RBPs/RPKMs/eCLIP-K562-all_RBPs_and_mocks-max_peaks-narrow-merged-coverage_intersect-narrow_clusters-ZRANB2.RPKM.tab","/media/skgthab/storage/UCL/2017.04.28-CLIPo-iONMF-ENCODE/CLIPo-eCLIP-K562-merged-narrow_clusters-RBPs_and_mock/RBPs/RPKMs/eCLIP-K562-all_RBPs_and_mocks-max_peaks-narrow-merged-coverage_intersect-narrow_clusters-all.tab")
'''
    
    