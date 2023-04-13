'''
Created on FEB 16 2016

@author: Nejc Haberman

Input should be in the folloqing order if you used this command: bedtools bamtobed -bedpe -mate1 -i data.bam > data.bed
And make sure that bam is sorted by name -n before running bedtools bamtobed

chr7    47565377        47565405        chr7    47565377        47565408        AATTT:SN1001:449:HGTN3ADXX:1:1110:20859:51900   255     +       -
chr14   102475740       102475779       chr14   102475729       102475773       AAATG:SN1001:449:HGTN3ADXX:1:1113:13799:66613   255     -       +
chr22   38881742        38881780        chr22   38881776        38881820        ATGTA:SN1001:449:HGTN3ADXX:2:1213:10755:25884   255     +       -
chrX    118553719       118553756       chrX    118553718       118553756       TTGCG:SN1001:449:HGTN3ADXX:1:1208:7358:45512    255     -       +
chr13   31038578        31038615        chr13   31038612        31038656        AAAGG:SN1001:449:HGTN3ADXX:1:2208:2880:92958    255     +       -
chr17   43012639        43012678        chr17   43012681        43012725        AAAAC:SN1001:449:HGTN3ADXX:1:1201:10544:64377   255     +       -
chr7    157731737       157731769       chr7    157731737       157731769       ACATC:SN1001:449:HGTN3ADXX:2:1209:17840:86837   255     +       -
chr7    138268719       138269499       chr7    138268718       138269499       TAGAG:SN1001:449:HGTN3ADXX:2:2213:19017:97248   255     -       +


'''


import sys

def BEDpaired2singleBED (fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        strand1 = col[8]
        strand2 = col[9]
        if strand2 == '+':
            start = col[4]
            end1 = col[2]
            end2 = col[5]
            fout.write(col[0] + '\t' + str(start) + '\t' + str(end1) + '\t' + str(end2) + '\t' + '' + '\t' + strand2 + '\n')
        elif strand2 == '-':
            start = col[5]
            end1 = col[1]
            end2 = col[4]
            fout.write(col[0] + '\t' + str(end1) + '\t' + str(start) + '\t' + str(end2) + '\t' + '' + '\t' + strand2 + '\n')
        line = fin.readline()
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDpaired2singleBED(fname_in, fname_out)
else:
    print("python BED2BEDgraph.py <input_file> <output_file>")

