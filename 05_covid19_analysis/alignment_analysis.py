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
import ast
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
    padding_novaseq = 17
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
    df_strains["Strain"] = [x.replace("-0","-") for x in df_strains["Strain"].values]

    if out_file != None:
        df_strains.to_csv(out_file,index=False)
    else:
        return df_strains

def tableGenomicMutationsMultiple(out_file):

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
    df.drop_duplicates("strain",inplace=True)
    print(df.shape)




    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2.drop(df2.index[0], inplace = True)

    ############ run gene wise
    # region_name ="S"
    # if region_name == "nsp12":
    #     region_start = 13442
    #     region_stop = 16237
    # elif region_name =="S":
    #     region_start = 21563
    #     region_stop = 25385
    

    # print("selecting region--"+region_name)
    # print(df2.shape)
    # df2 = df2.iloc[region_start:region_stop,:]
    # print(df2.shape)
    ### to calculate muations per strain in genome


    mismatch =[]
    strain_counter = 0
    mismatch.append(["Reference",0,[]])
    for column in df2.columns[1:]:
        mutation_counter = 0
        mutation = []
        for indx,ref,alt in zip(df2.index.values,df2["MN908947"].values,df2[column].values):
            if ref != alt and alt in ['a','t','g','c']:
                mutation.append(ref.upper()+str(indx)+alt.upper())
                mutation_counter += 1
        mismatch.append([column,mutation_counter,mutation])
        
        strain_counter += 1
        if strain_counter % 1000 == 0:
            print("genomewide mutation table completed for.."+str(strain_counter)+" strains")

    df_strains = pd.DataFrame(mismatch)
    df_strains.columns = ["Strain","mutation_count","mutation_info"]
    df_strains["Strain"] = [x.replace("-0","-") for x in df_strains["Strain"].values]
    # df_strains.to_csv("muttable.csv",index=False)  

    if out_file != None:
        df_strains.to_csv(out_file,index=False)
    else:
        return df_strains


def queryRegion(alignment_file,out_file):

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

    ## loof for region of interest    
    # sample_seq =  "".join (s for s in df.iloc[0,212:])
    # sample_seq[27628:27994]
    
    ## confirm 
    # df.iloc[0,27840:28206]

    df_roi = df.iloc[:,27840:28206]

    df_roi = df.iloc[:,27840:28230]


    df2 = df_roi.apply(counvals,axis=1)
    df2["Strain"] = df[0]

    df2.to_csv("orf8_v2.tsv",sep='\t')

def counvals(df):
    df['n']=np.sum(df=='n')
    df['-']=np.sum(df=='-')
    return df

def AnnotateTabledGenomicMutations(genomic_mutation_table,out_file):

    gene_name={
    'S':'21563-25384',
    'E':'26245-26472',
    'M':'26523-27191',
    'N':'28274-29533',
    'NSP12':'13446-16236'}

    df = pd.read_csv(genomic_mutation_table)

    mutation_wise =[]
    for indx,row in df.iterrows():
        for mutation in ast.literal_eval(row["mutation_info"]):
            position = int(mutation.replace("A","").replace("T","").replace("G","").replace("C",""))
            gene = ""

            for g in gene_name:
                start = int(gene_name[g].split('-')[0])
                end = int(gene_name[g].split('-')[1])

                if position>=start and position<=end:
                    gene = g
            mutation_wise.append([row["Strain"],row["mutation_count"],mutation,gene, row["mutation_info"]])

    df_mutation_wise = pd.DataFrame(mutation_wise)
    df_mutation_wise.columns = ["Strain","mutation_count","mutation","gene","all_mutations"]
    df_mutation_wise.to_excel(out_file,index=False)

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

def getStrainswithPosition(alignment_file,region,position,novaseq_adj=None,out_file=None):
    region_start = region[0]
    region_stop = region[1]

    filt = []
    #### read alignment files for reference
    alignment = AlignIO.read("data/4_15/Houston.4-14.clean.fa","fasta")
    # alignment = AlignIO.read("../../data/4_15/Houston.4-14.clean.fa","fasta")

    for line in alignment:
        temp =[]
        temp.append(line.id)
        print(line.id)
        for nt in str(line.seq[region_start:region_stop]):
            temp.append(nt.lower())
        filt.append(temp)
        break

    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    ## adjustment for novaseq
    if novaseq_adj == None:
       novaseq_adj = 0


    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq[region_start+novaseq_adj:region_stop+novaseq_adj]):
            temp.append(nt.lower())
        filt.append(temp)

    df = pd.DataFrame(filt)
    # df[0]=[x.replace("V-",".").replace("-","").replace(".","-") for x in df[0].values]
    df.drop_duplicates(0,inplace=True)

    print("before removing ...")
    print(df.shape)

    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

    region_start = 21563
    region_stop = 25385

    col = []
    col.append("id")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col

    return df.ix[:,['id',23403]]

def tableGeneMutations(alignment_file,region,region_name,out_file):

    region_start = region[0]
    region_stop = region[1]

    filt = []
    #### read alignment files for reference
    alignment = AlignIO.read("data/4_15/Houston.4-14.clean.fa","fasta")
    # alignment = AlignIO.read("../../data/4_15/Houston.4-14.clean.fa","fasta")

    for line in alignment:
        temp =[]
        temp.append(line.id)
        print(line.id)
        for nt in str(line.seq[region_start:region_stop]):
            temp.append(nt)
        filt.append(temp)
        break

    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []

    ## adjustment for novaseq
    novaseq_adj = 17
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq[region_start+novaseq_adj:region_stop+novaseq_adj]):
            temp.append(nt)
        filt.append(temp)

    df = pd.DataFrame(filt)
    print(df.shape)
    print(df.head())


    if region_name == "nsp12":
        region_start = 13442
        region_stop = 16237
    elif region_name =="S":
        region_start = 21563
        region_stop = 25385

    
    col = []
    col.append("strain")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col

    print(df.shape)
    print(df.head())

    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]

    print(df2.shape)
    print(df2.head())
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
    df_strains["Strain"] = [x.replace("-0","-") for x in df_strains["Strain"].values]

    if out_file != None:
        df_strains.to_csv(out_file,index=False)
    else:
        return df_strains


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
    else:
        region=[21508,25330]
        position=23403

        alignment_file ="/home/tmhsxs240/COVID_19/data/8_11_all_5085/Houston.Aug12.clean.fa"
        df_5805 = getStrains(alignment_file,out_file=None)
        df_5805["group"]="g1_5805"

        alignment_file ="/home/tmhsxs240/COVID_19/data/11_20/Nov-19.trimmed.clean.fa"
        df_oct = getStrains(alignment_file,out_file=None)
        df_oct["group"]="g2_NOVA1_4"

        alignment_file ="/home/tmhsxs240/COVID_19/data/12_1/Nov-29.trimmed.clean.fa"
        df_nova1 = getStrains(alignment_file,out_file=None)
        df_nova1["group"]="g3_NOVA_R5"


        df_all = df_5805.append(df_oct)
        df_all = df_all.append(df_nova1)


            ## remove 1255 and 1343 and two training runs
        df_all = df_all[df_all["id"] != "MCoV-1255"]
        df_all = df_all[df_all["id"] != "MCoV-1343"]
        df_all = df_all[df_all["id"] != "MCoV-743"]
        df_all = df_all[df_all["id"] != "MCoV-1219"]

        ### more noise data from 5x depth
        df_all = df_all[df_all["id"] != "MCoV-12783"]
        df_all = df_all[df_all["id"] != "MCoV-12765"]
        df_all = df_all[df_all["id"] != "MCoV-12805"]

        df_all.sort_values("group",inplace=True)
        df_all.drop_duplicates("id",inplace=True)
        df_all.columns = ["Strain","Variant614","Group"]
        df_all["Variant614"] = ["D614" if x=="a" else "G614" if x=="g" else "N" for x in df_all["Variant614"].values ]
        df_all["Strain"] = [x.replace("-0","-") for x in df_all["Strain"].values]    
        
        ref_dir="/home/tmhsxs240/COVID_19/reference/"
        db= pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains.xlsx")

        dfjoin = pd.merge(df_all,db,how="left",on="Strain",indicator=True)
        dfjoin = dfjoin[dfjoin["_merge"]=="both"]
        dfjoin.to_csv("data/9_29_allD614G/D614G_good_strains_from_alignment_until_nova_r3.csv",index=False)


