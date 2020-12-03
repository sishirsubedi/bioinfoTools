import pandas as pd
import os
import re
import ast
from Bio import AlignIO

def combineGenome():
    print("processing..5085")
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

    alignment_file ="/home/tmhsxs240/COVID_19/data/8_11_all_5085/Houston.Aug12.clean.fa"
    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt.lower())
        filt.append(temp)

    df = pd.DataFrame(filt)
    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

    region_start = 266
    region_stop = 266+(len(df.columns)-1)

    col = []
    col.append("strain")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col

    print(df.shape)
    print(df.head())

    print("processing previous nta")
    filt = []
    padding_novaseq = 16
    #### read alignment files for reference
    alignment = AlignIO.read("/home/tmhsxs240/COVID_19/data/4_15/Houston.4-14.clean.fa","fasta")
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


    alignment_file ="/home/tmhsxs240/COVID_19/data/11_20/Nov-19.trimmed.clean.fa"
    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt.lower())
        filt.append(temp)


    df_p2 = pd.DataFrame(filt)
    print(df_p2.shape)
    print(df_p2.head())

    ### remove any "-" from the reference
    df_p2 = df_p2.drop(df_p2.columns[df_p2.iloc[0,:] =='-'],axis=1)
    df_p2 = df_p2.drop(df_p2.columns[df_p2.iloc[0,:].isnull()],axis=1)

    region_start = 266
    region_stop = 266+(len(df_p2.columns)-1)

    col = []
    col.append("strain")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df_p2.columns = col

    print(df_p2.shape)
    print(df_p2.head())

    print("processing new nta")
    filt = []
    padding_novaseq = 26
    #### read alignment files for reference
    alignment = AlignIO.read("/home/tmhsxs240/COVID_19/data/4_15/Houston.4-14.clean.fa","fasta")
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


    alignment_file ="/home/tmhsxs240/COVID_19/data/12_1/Nov-29.trimmed.clean.fa"
    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt.lower())
        filt.append(temp)


    df_p3 = pd.DataFrame(filt)
    print(df_p3.shape)
    print(df_p3.head())

    ### remove any "-" from the reference
    df_p3 = df_p3.drop(df_p3.columns[df_p3.iloc[0,:] =='-'],axis=1)
    df_p3 = df_p3.drop(df_p3.columns[df_p3.iloc[0,:].isnull()],axis=1)

    region_start = 266
    region_stop = 266+(len(df_p3.columns)-1)

    col = []
    col.append("strain")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df_p3.columns = col

    print(df_p3.shape)
    print(df_p3.head()) 


    ################combine all #################

    print(df.shape)
    df = df.append(df_p2)
    df.reset_index(inplace=True,drop=True)

    print(df.shape)
    df = df.append(df_p3)
    df.reset_index(inplace=True,drop=True)

    print(df.shape)
    df["strain"] = [x.replace("-0","-") for x in df["strain"].values]
    df.drop_duplicates("strain",inplace=True)
    print(df.shape)

    return df


#### common snps
def run():
    
    df = combineGenome()

    dfselected = df[["strain",23403,21707,24718,24026,21575,23481,23604,21638,23224,24812]]

    dfselected.rename(columns={
    23403:"D614G",
    21707:"H49Y",
    24718:"F1052L",
    24026:"L822F",
    21575:"L5F",
    23481:"S640F",
    23604:"P681H",
    21638:"P26S",
    23224:"E554D",
    24812:"D1084Y"
    },inplace=True)

    dfselected["D614G"] =["D" if x =="a" else  "G" if x =="g" else "OTHER" for x in dfselected["D614G"].values]
    dfselected["H49Y"] =["H" if x =="c" else  "Y" if x =="t" else "OTHER" for x in dfselected["H49Y"].values]
    dfselected["F1052L"] =["F" if x =="c" else  "L" if x =="a" else "OTHER" for x in dfselected["F1052L"].values]
    dfselected["L822F"] =["L" if x =="c" else  "F" if x =="t" else "OTHER" for x in dfselected["L822F"].values]
    dfselected["L5F"] =["L" if x =="c" else  "F" if x =="t" else "OTHER" for x in dfselected["L5F"].values]
    dfselected["S640F"] =["S" if x =="c" else  "F" if x =="t" else "OTHER" for x in dfselected["S640F"].values]
    dfselected["P681H"] =["P" if x =="c" else  "H" if x =="a" else "OTHER" for x in dfselected["P681H"].values]
    dfselected["P26S"] =["P" if x =="c" else  "S" if x =="t" else "OTHER" for x in dfselected["P26S"].values]
    dfselected["E554D"] =["E" if x =="g" else  "D" if x =="t" else "OTHER" for x in dfselected["E554D"].values]
    dfselected["D1084Y"] =["D" if x =="g" else  "Y" if x =="t" else "OTHER" for x in dfselected["D1084Y"].values]


    ## remove 1255 and 1343 and two training runs
    # "MCoV-1255","MCoV-1343","MCoV-743","MCoV-1219"
    ### more noise data from 5x depth
    # "MCoV-12783","MCoV-12765","MCoV-12805"
    ###more training runs not present in soft database
    # "MCoV-12123","MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867",
    ###not found in soft
    # "MCoV-11867"



    remove_index = df[df["Strain"].isin(["MCoV-1255","MCoV-1343","MCoV-743","MCoV-1219",
    "MCoV-12783","MCoV-12765","MCoV-12805",
    "MCoV-12123","MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867","MCoV-11867"])].index
    df.drop(remove_index,inplace=True)

    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    dflog = pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains.xlsx")

    dfselected.rename(columns={'strain':'Strain'},inplace=True)
    dfjoin = pd.merge(dfselected,dflog,on="Strain",how="left",indicator=True)
    dfjoin = dfjoin[dfjoin._merge=="both"]
    dfjoin.to_excel(ref_dir+"All_Strains_common_SNPs.xlsx",index=False)