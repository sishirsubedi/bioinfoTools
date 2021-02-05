import pandas as pd
import os
from Bio import AlignIO
import gen_utils.gen_covid as covid


def genomeSequenceFromAlignment(alignment_file):
    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        strain_info =[]
        strain_info.append(strain.id)
        for nt in str(strain.seq):
            strain_info.append(nt.lower())
        sequences.append(strain_info)
    return sequences

def genomeSequenceFromAlignmentwithReferencePadding(reference_alignment,alignment_file,out_file):

    filt = []
    padding_novaseq = 17
    #### read alignment files for reference
    alignment = AlignIO.read(reference_alignment,"fasta")
    for line in alignment:
        temp =[]
        temp.append(line.id)
        print(line.id)
        for i in range(padding_novaseq):
            temp.append("-")
        for nt in str(line.seq):
            temp.append(nt)

        filt.append(temp)
        break

    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt)
        filt.append(temp)


    return filt

def geneSequenceFromAlignment(alignment_file,gene,coordinates,reference=False,adjustment=0):

    gene_start = int(coordinates.split(':')[0])
    gene_stop = int(coordinates.split(':')[1])

    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        strain_info =[]
        strain_info.append(strain.id)
        for nt in str(strain.seq[gene_start+adjustment:gene_stop+adjustment]):
            strain_info.append(nt.lower())
        sequences.append(strain_info)
        if reference:
            break

    return sequences

def strainFromAlignment(alignment_file,runid):

    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        sequences.append([strain.id,runid])

    return sequences


def primerSequencesFromAlignment(primer_pairs,alignment_file):

    # primer_pairs = [[15513:15535],[15551:15577],[15587:15613]]

    all_alignment = AlignIO.read(alignment_file,"fasta")

    all = []
    for line in all_alignment:
        temp =SeqRecord(Seq(str(line.seq[primer_pairs[0]])+"-"\
                            +str(line.seq[primer_pairs[1]])+"-"\
                            +str(line.seq[primer_pairs[2]])),line.id,line.name,line.description)
        all.append(temp)
    alignment = MultipleSeqAlignment(all)

    filt = []
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt)
        filt.append(temp)

    return filt


def analyzeGenomicRegionwithReference(reference_alignment,alignment_file,region_name):

    if region_name == "S":
        region="21508:25330"
    elif region_name == "nsp12":
        region="13387:16182"

    reference_strain = geneSequenceFromAlignment(reference_alignment,region_name,region,reference=True)
    reference_sequence = "".join(x for x in reference_strain[0][1:])

    print("reference sequence for "+region_name)
    print(reference_sequence[0:100])
    print(reference_sequence[-100:])


    all_strains = genomeSequenceFromAlignment(alignment_file)
    df = pd.DataFrame(all_strains)

    query_sequence = "".join(x for x in df[df[0]=='NC_045512.2'].values[0][1:])

    start = query_sequence.find(reference_sequence[0:30])+1
    stop = start + len(reference_sequence)

    ##adjustment for reference just in case there are "gaps (-)" inserted in the reference during alignment
    ## we need to remove those gaps and adjust genomic coordinates
    adjustments_for_reference = collections.Counter(df[df[0]=='NC_045512.2'].values[0][start:stop])['-']
    stop += adjustments_for_reference

    sel_columns = [0]
    for i in range(start,stop):sel_columns.append(i)
    df = df.iloc[:,sel_columns]

    df.sort_index(ascending=False,inplace=True)
    df.reset_index(drop=True,inplace=True)
    df.iloc[0,0] = covid.covid_basics.reference

    return df


def get_df_from_multiple_alignments(all_aligned_fasta_files,out_dir):

    df_list = []
    for alignment_file in all_aligned_fasta_files:

        all_strains = genomeSequenceFromAlignment(out_dir+alignment_file)
        df = pd.DataFrame(all_strains)

        df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

        sel_columns = []
        for i in range(0,29904):sel_columns.append(i)
        df.columns = sel_columns
        df.loc[0,0] = covid.covid_basics.reference

        df_list.append(df)

    df = df.append(df_list[0])
    for x in range(1,len(df_list)):
        df = df.append(df_list[x])

    df.drop_duplicates(0,inplace=True)
    df.reset_index(inplace=True,drop=True)

    return df
