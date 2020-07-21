""" analyze alignment file

##### example:

 /opt/python3/bin/python3 alignment_analysis.py \
 --request "get-strains" \
 --alignment "/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa" \
 --output "/home/tmhsxs240/COVID_19/data/Ns_Strain_Comparison/D614G_Npercent_analysis.csv"

"""

import pandas as pd
import numpy as np
import argparse
import os
import re
from Bio import AlignIO, SeqIO, Entrez

def getStrains(alignment_file,out_file=None):
    region_start = 0
    region_stop = 1
    filt = []
    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq[region_start:region_stop]):
            temp.append(nt)
        filt.append(temp)
    df_combine = pd.DataFrame(filt)
    df_combine.columns = ["Strain","Seq"]

    df_combine["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("v3","").replace("V3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("_","-") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-0","-") for x in df_combine["Strain"].values]
    df_combine.drop_duplicates("Strain",inplace=True)

    if out_file != None:
        df_combine.to_csv(out_file,index=False)
    else:
        return df_combine

def tableGenomicMutations(alignment_file,out_file):

    filt = []
    #### read alignment files for reference
    alignment = AlignIO.read("/home/tmhsxs240/COVID_19/data/4_15/Houston.4-14.clean.fa","fasta")
    for line in alignment:
        temp =[]
        temp.append(line.id)
        print(line.id)
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

    df = pd.DataFrame(filt)
    print(df.shape)
    print(df.head())

    ###keep inhouse samples only, remove GISAID
    df = df[~df[0].str.contains("EPI_ISL")]

    ### remove any "-" from the reference
    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

    region_start = 266
    region_stop = 266+(len(df.columns)-1)

    col = []
    col.append("strain")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col

    print(df.shape)
    print(df.head())

    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]

    ### to calculate muations per strain in genome

    mismatch =[]
    for column in df2:
        if column =="MN908947":
            mismatch.append([column,0,[]])
        else:
            counter = 0
            mutation = []
            for indx,row in df2[column].items():
                if df2.ix[indx,column]=='a' or df2.ix[indx,column]=='t' or df2.ix[indx,column]=='g' or df2.ix[indx,column]=='c':
                    if df2.ix[indx,column] != df2.ix[indx,"MN908947"] :
                        change = df2.ix[indx,"MN908947"].upper()+str(indx)+df2.ix[indx,column].upper()
                        mutation.append(change)
                        counter += 1
            mismatch.append([column,counter,mutation])

    df_strains = pd.DataFrame(mismatch)
    df_strains.columns = ["Strain","mutation_count","mutation_info"]

    if out_file != None:
        df_strains.to_csv(out_file,index=False)
    else:
        return df_strains


def primerAlignment(alignment_file,out_file):
    #### read alignment files
    all_alignment = AlignIO.read(alignment_file,"fasta")

    all = []
    for line in all_alignment:
        temp =SeqRecord(Seq(str(line.seq[15513:15535])+"-"+str(line.seq[15551:15577])+"-"+str(line.seq[15587:15613])),line.id,line.name,line.description)
        all.append(temp)
    alignment = MultipleSeqAlignment(all)

    filt = []
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt)
        filt.append(temp)
    df = pd.DataFrame(filt)
    mismatch =[]
    for column in df:
        df_clean = df[df[column]!='n']
        df_clean = df_clean[df_clean[column]!='k']
        mismatch.append(sum(df_clean[column]!= df_clean.iloc[0,column]))

    col = []
    col.append("id")
    for x in range(1,len(mismatch),1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col
    df = df.T
    df.rename(columns=df.iloc[0],inplace=True)
    df=df.iloc[1:,:]
    df["mismatch_count"] = mismatch[1:]
    df = df[df["mismatch_count"]>=1]
    sample_number = 5460
    variants =[]
    variant_count =[]
    for indx,row in df.iterrows():

        vc = str(row[0:sample_number].value_counts()).split("\n")
        variant_count.append(vc[:len(vc)-1])

        t = row[0:sample_number].unique()
        variants.append([x for x in t if x!= row[0]])

    df["variants"] = variants
    df["variant_count"] = variant_count
    df.reset_index(inplace=True)
    df = df[['index', 'MN908947.3','mismatch_count', 'variants','variant_count']]

    if out_file != None:
        df.to_excel(out_file,index=False)
    else:
        return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign mutation status, sequence quality, and alignment status to strain",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--request", help="task to perform")
    parser.add_argument("--alignment", help="alignment file")
    parser.add_argument("--output", help="output file")
    args = parser.parse_args()

    if args.request == "get-strains":
        getStrains(args.alignment,args.output)
    elif args.request == "mutations-table":
        tableGenomicMutations(args.alignment,args.output)
