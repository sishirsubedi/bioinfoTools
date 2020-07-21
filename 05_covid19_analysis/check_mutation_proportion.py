import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import seaborn as sns
from Bio import AlignIO, SeqIO, Entrez
import numpy as np
import os
import re
import ast
import alignment_analysis

alignment_file =""
genomic_mutation_table = "1_StrainsAndMutations_entire_genome.csv"
alignment_analysis.tableGenomicMutations(alignment_file,genomic_mutation_table):


df = pd.read_csv(genomic_mutation_table)
df["Strain"] = [x.replace("-0","-") for x in df["Strain"].values]
df = df[df.Strain!="MN908947"]
df = df[df.Strain!="MCoV-1343"]
df = df[df.Strain!="MCoV-1255"]


# df_sample_log= pd.read_excel("3_MCoV Sample Log 6-22-20 for Paul.xlsx")

df_sample_log= pd.read_excel("MCoV Metadata for PATRIC 6-23.xlsx")
df_sample_log = df_sample_log[['MRN', 'Full Order Number', 'Musser Lab No.', 'GENDER','AGE','DeceasedYN-June22','ETHNIC_GROUP']]
df_sample_log.columns = ['MRN', 'ORDER_ID', 'MCoVNumber','GENDER','AGE','DEATH','ETHNIC_GROUP']


df_db = pd.read_csv("1_Curated_Covid_database.csv")
df_db.COLLECTION_DT = pd.to_datetime(df_db.COLLECTION_DT)
df_db.VERIFIED_DT = pd.to_datetime(df_db.VERIFIED_DT)

dfjoin_mrn = pd.merge(df_sample_log,df_db, on="ORDER_ID",how="left",indicator=True)
dfjoin_mrn = dfjoin_mrn[dfjoin_mrn._merge=="both"]
dfjoin_mrn = dfjoin_mrn[['MRN_x','ORDER_ID', 'MCoVNumber', 'COLLECTION_DT','VERIFIED_DT','GENDER','AGE','DEATH','ETHNIC_GROUP']]
dfjoin_mrn.columns = ['MRN','ORDER_ID', 'MCoVNumber', 'COLLECTION_DT','VERIFIED_DT','GENDER','AGE','DEATH','ETHNIC_GROUP']

# dfjoin_mrn = dfjoin_mrn[dfjoin_mrn["MCoVNumber"].notnull()]
# dfjoin_mrn.sort_values("COLLECTION_DT",inplace=True)
# dfjoin_mrn.drop_duplicates(["MCoVNumber"],inplace=True)

method="COLLECTION_DT"
dffinal = pd.merge(df,dfjoin_mrn, left_on="Strain",right_on="MCoVNumber",how="left",indicator=True)
dffinal = dffinal[dffinal._merge=="both"]
dffinal.set_index(method,inplace=True)

wave=[]
for indx,row in dffinal.iterrows():
    if indx < pd.Timestamp('2020-05-12 00:00:00'):
        wave.append(1)
    else:
        wave.append(2)

dffinal['wave'] =wave

df_week = dffinal.to_period(freq='W-MON')
df_week.reset_index(inplace=True)

moi = ["A23403G","C24034T","C23191T","T24982C","C24718A"]
for mutation in moi:
    mut_614=[]
    for indx,row in df_week.iterrows():
        mut_list = ast.literal_eval(row.mutation_info)
        if mutation in mut_list:
            mut_614.append(1)
        else:
            mut_614.append(0)
    df_week[mutation] = mut_614

###group age and gender


df_week.to_csv("df_week.csv")
# df_week.D614G.value_counts()


df_week.groupby(method).count()
df_week.groupby(method).agg(sum)

df_week.groupby(method)["AGE"].mean()
df_week.groupby(method)["GENDER"].value_counts()
df_week.groupby(method)["DEATH"].value_counts()


df_week[df_week.wave==1].AGE.plot.hist()
plt.title("Age distribution - wave 1")
plt.xlabel("Age")
plt.savefig("age_dist_w1.png");plt.close()

df_week[df_week.wave==2].AGE.plot.hist()
plt.title("Age distribution - wave 2")
plt.xlabel("Age")
plt.savefig("age_dist_w2.png");plt.close()

###############
