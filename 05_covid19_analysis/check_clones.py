import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO, Entrez
import ast
import pylab
import pandas as pd
import seaborn as sns
import numpy as np
import os
import re
import alignment_analysis

alignment_analysis.tableGenomicMutations("/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa","5_StrainsAndMutations_entire_genome.csv")
df_mutations_per_strain = pd.read_csv("5_StrainsAndMutations_entire_genome.csv")



df= pd.read_excel("COVID_SNP_MAP_S_NTA_variants_7_2_updated.xlsx")
# df = df[df["NTA_Mismatch_Count"]>1]
# df = df[df["NTA_Mismatch_Count"]<10]


df = df[['NTA_Genomic_Locus', 'NTA_Referene', 'NTA_Mismatch_Count', 'NTA_Gene','NTA_Variants', 'NTA_Variants_Info']]

combine =[]
for indx,row in df.iterrows():
    alt_variant = ast.literal_eval(row["NTA_Variants"])[0]
    variants = ast.literal_eval(row["NTA_Variants_Info"])
    for variant in variants:
        if variant == alt_variant:
            combine.append([row["NTA_Mismatch_Count"],row["NTA_Referene"].upper()+str(row['NTA_Genomic_Locus'])+alt_variant.upper(),variants[variant][1]])
df_combine = pd.DataFrame(combine)
df_combine.columns =[ 'NTA_Mismatch_Count','NTA_Variants', 'NTA_Variants_Info']

df_combine["count_len"] = [len(x) for x in df_combine["NTA_Variants_Info"]]

df_combine  = df_combine[df_combine["count_len"]==df_combine["NTA_Mismatch_Count"]]


df_combine.NTA_Mismatch_Count.value_counts()
df_combine = df_combine[df_combine["NTA_Mismatch_Count"]>1]

clone_match={}
for num in df_combine.NTA_Mismatch_Count.unique():

    clone_match[num]=0
    df_selected = df_combine[df_combine["NTA_Mismatch_Count"]==num]

    all_strains = []
    for indx,strain_list in df_selected.NTA_Variants_Info.items():
        for strain in strain_list:
            if strain not in all_strains:
                all_strains.append(strain)

    repeat_strains = {}
    for s in all_strains: repeat_strains[s] =0

    for indx,strain_list2 in df_selected.NTA_Variants_Info.items():
        for strain2 in all_strains:
            if strain2 in strain_list2:
                repeat_strains[strain2] += 1

    for s in repeat_strains:
        if repeat_strains[s] >1:
            clone_match[num] += 1


df_combine.to_excel("Clonal_analysis_from_spike_protein.xlsx",index=False)

df_db = pd.read_csv("1_StrainsAndMutations_entire_genome.csv")
clones = df_combine[df_combine.NTA_Variants=="C22444T"]["NTA_Variants_Info"].values[0]
df_db = df_db[df_db.Strain.isin(clones)]
df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_curated = df_curated[['MCoVNumber', 'ORDER_ID', 'MRN']]
df_curated.columns =['Strain', 'ORDER_ID','MRN']
dfjoin_selected = pd.merge(df_db,df_curated,how="left",on="Strain",indicator=True)
dfjoin_selected.to_excel("7_isolates.xlsx",index=False)

df_db = pd.read_csv("1_StrainsAndMutations_entire_genome.csv")
clones = df_combine[df_combine.NTA_Variants=="C21707T"]["NTA_Variants_Info"].values[0]
df_db = df_db[df_db.Strain.isin(clones)]
df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_curated = df_curated[['MCoVNumber', 'ORDER_ID', 'MRN']]
df_curated.columns =['Strain', 'ORDER_ID','MRN']
dfjoin_selected = pd.merge(df_db,df_curated,how="left",on="Strain",indicator=True)
dfjoin_selected.to_excel("3_A_isolates.xlsx",index=False)

df_db = pd.read_csv("1_StrainsAndMutations_entire_genome.csv")
clones = df_combine[df_combine.NTA_Variants=="A23083G"]["NTA_Variants_Info"].values[0]
df_db = df_db[df_db.Strain.isin(clones)]
df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_curated = df_curated[['MCoVNumber', 'ORDER_ID', 'MRN']]
df_curated.columns =['Strain', 'ORDER_ID','MRN']
dfjoin_selected = pd.merge(df_db,df_curated,how="left",on="Strain",indicator=True)
dfjoin_selected.to_excel("3_B_isolates.xlsx",index=False)

#################### add zip code information

df_db = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_zipcode = pd.read_csv("mrn_to_zipcode.csv")
df_zipcode.ZIP = [str(x)[0:5] for x in df_zipcode.ZIP]
df_db_join = pd.merge(df_db,df_zipcode,on="MRN",how="left",indicator=True)
df_db_join = df_db_join[['MCoVNumber', 'ORDER_ID', 'MRN', 'COLLECTION_DT', 'ZIP']]
df_db_join.columns = ['Strain', 'ORDER_ID', 'MRN', 'COLLECTION_DT', 'ZIP']
df_db_join.to_excel("4_Curated_MCOV_MRN_Strains_with_zipcode.xlsx")

df_alignment = alignment_analysis.getStrains("/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa")

dfjoin_selected = pd.merge(df_alignment,df_db_join,how="left",on="Strain",indicator=True)
dfjoin_selected = dfjoin_selected[dfjoin_selected._merge=="both"]
dfjoin_selected.drop_duplicates("Strain",inplace=True)


df_clonal = pd.read_excel("Clonal_analysis_from_spike_protein.xlsx")

zip_codes =[]
mrn_codes =[]
for indx,row in df_clonal.iterrows():
    zip = []
    mrn = []
    strain_list = ast.literal_eval(row.NTA_Variants_Info)
    for strain in strain_list:
        print(strain)
        strain = strain.replace("-0","-")
        if dfjoin_selected[dfjoin_selected.Strain==strain].shape[0] != 0:
            zip.append(dfjoin_selected[dfjoin_selected.Strain==strain]["ZIP"].values[0])
            mrn.append(dfjoin_selected[dfjoin_selected.Strain==strain]["MRN"].values[0])
        else:
            zip.append(0)
            mrn.append(0)
    zip_codes.append(zip)
    mrn_codes.append(mrn)

df_clonal["zip_codes"] = zip_codes

df_clonal["MRNs"] = mrn_codes

df_clonal["zip_codes_counts"] = [len(set(x)) for x in df_clonal["zip_codes"]]

df_clonal["MRN_counts"] = [len(set(x)) for x in df_clonal["MRNs"]]

df_clonal.to_excel("Clonal_analysis_from_spike_protein_with_zipcode.xlsx")


#############houston border #################
############################################

BBox = (-96.3940,-94.5813,30.3551,29.3247)

from matplotlib.backends.backend_pdf import PdfPages
from uszipcode import SearchEngine

#####plot zipcodes

df_clonal_mini = df_clonal[df_clonal.NTA_Mismatch_Count>=0]

df_clonal_mini.sort_values("NTA_Mismatch_Count",ascending=False,inplace=True)


search = SearchEngine(simple_zipcode=True)


with PdfPages('multipage_pdf_all.pdf') as pdf:
    for indx,row in df_clonal_mini.iterrows():

        fig, ax = plt.subplots(figsize = (15,10))
        for zcode in row.zip_codes:

            # zcdb = ZipCodeDatabase()

            if zcode is None or zcode=='nan': continue
            print(zcode)
            try:
                zipcode = search.by_zipcode(int(zcode))

                ax.scatter(zipcode.lng, zipcode.lat, zorder=1, c='r', s=50)
            except:
                print("failed---"+str(zcode))


        ax.set_title(row.NTA_Variants+" , total isolates- "+str(row.NTA_Mismatch_Count)+" , total patient- "+str(row.MRN_counts)+ " ,total zipcodes- "+str(row.zip_codes_counts))
        ax.imshow(ruh_m, zorder=0, extent = BBox, aspect= 'equal')
        pdf.savefig(fig)
