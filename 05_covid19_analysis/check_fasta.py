
import pandas as pd
from Bio import AlignIO, SeqIO, Entrez
import re
import sys

def getsubstring(position,mode):
    if mode=='query':
        ref_seq = []
        for feature in SeqIO.parse("reference/nCoV-2019.reference.fasta", "fasta"):
            ref_seq = feature.seq
        return ref_seq[position-20:position-1]
    if mode=='seq':
        ref_seq = []
        for feature in SeqIO.parse("reference/nCoV-2019.reference.fasta", "fasta"):
            ref_seq = feature.seq
        return ref_seq[position-20:position+10]


##################################################
# sample =""
sample=sys.argv[1]
position=int(sys.argv[2])

#############################################################
query_string = getsubstring(position-1,'query')

sample_seq = []
for feature in SeqIO.parse(sample, "fasta"):
    sample_seq = feature.seq

sample_location = re.search(str(query_string), str(sample_seq)).span()

print('     ')
print('     ')
print("Sample    : ",str(sample.split('/')[3]))
print("Position    : ",str(position))
print("Reference query    : ",str(query_string))
print("Sample query match : ",sample_seq[sample_location[0]:sample_location[1]])

print('                            :  --------------------v----------')
print("Reference seq               : ",str(getsubstring(position-1,'seq')))
print("sample seq with mutation    : ",sample_seq[sample_location[0]:sample_location[1]+10])
print('     ')
print('     ')
