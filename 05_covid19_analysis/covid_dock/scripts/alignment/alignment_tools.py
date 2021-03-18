import pandas as pd
import os
from Bio import AlignIO
import collections
import gen_utils.gen_covid as covid

def strainFromAlignment(alignment_file,runid):

    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        sequences.append([strain.id,runid])

    return sequences

def genomicSequenceFromAlignment(alignment_file,reference=False,adjustment=0):

    print("processing: "+alignment_file)

    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        strain_info =[]
        strain_info.append(strain.id)
        for nt in str(strain.seq[adjustment:]):
            strain_info.append(nt.lower())
        sequences.append(strain_info)
        if reference:
            break

    return sequences

def analyzeEntireGenomewithReference(alignment_file):

    all_strains = genomicSequenceFromAlignment(alignment_file)
    df = pd.DataFrame(all_strains)

    df.sort_index(ascending=False,inplace=True)
    df.reset_index(drop=True,inplace=True)
    df.iloc[0,0] = covid.covid_basics.reference

    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

    sel_columns = []
    for i in range(0,29904):sel_columns.append(i)
    df.columns = sel_columns

    return df

