import pandas as pd
from Bio import AlignIO, SeqIO, Entrez
import re
import sys
import ast
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



########################################################

ref_dir="/home/tmhsxs240/COVID_19/reference/"
print("running....genome wide mutation table")
alignment_analysis.tableGenomicMutations("/home/tmhsxs240/COVID_19/data/8_11/Houston.Aug12.clean.fa",ref_dir+"5_StrainsAndMutations_entire_genome.csv")
print("completed....")


###############################################################################



import pandas as pd
import vcf_analysis
import fasta_analysis

### get all data 
cd /home/tmhsxs240/COVID_19/data/9_29_allD614G/all_vcfs
ln -s ../../8_11_all_5085/8_11_all_vcfs/* .    ## 5085
ln -s ../../9_18_all_july/9_18_vcfs/* .        ## july 80-109
ln -s ../../10_23_all_augsept/10_23_vcfs/* .        ## sept 110-119

vcfs_directory = "data/9_29_allD614G/all_vcfs/"
position = 23403
out_file = "data/9_29_allD614G/D614G_all_isolates.csv"

vcf_analysis.analyzePositionVariant(vcfs_directory,position,out_file)
df = pd.read_csv("data/9_29_allD614G/D614G_all_isolates.csv")
df["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df["Strain"].values]
df["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df["Strain"].values]
df["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df["Strain"].values]
df["Strain"] = [x.replace("v3","").replace("V3","") for x in df["Strain"].values]
df["Strain"] = [x.replace("_","-") for x in df["Strain"].values]
df["Strain"] = [x.replace("-0","-") for x in df["Strain"].values]
df["Strain"] = [x.replace("-R","") for x in df["Strain"].values]
df.drop_duplicates("Strain",inplace=True)


df = df[df.Strain!="MCoV-12291"]
df = df[df.Strain!="MCoV-11867"]
df = df[df.Strain!="MCoV-743"]
df = df[df.Strain!="MCoV-1219"]


ref_dir="/home/tmhsxs240/COVID_19/reference/"
db= pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains.xlsx")

dfjoin = pd.merge(df,db,how="left",on="Strain",indicator=True)
dfjoin.to_csv("data/9_29_allD614G/D614G_all_withmeta.csv",index=False)



### get all data 
cd /home/tmhsxs240/COVID_19/data/9_29_allD614G/all_fastas
cp -avr /net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/Ns_Strain_Comparison/6_24_fastas/* .
cp -avr /net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/7_15/7_15_fastas/*  .
cp -avr /net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/8_11_all_5085/8_11_fastas/*  .
cp -avr /net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/9_18_all_july/9_18_fastas/*  .
cp -avr /net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/10_23_all_augsept/10_23_fastas/*  .


df_vcf = pd.read_csv("data/9_29_allD614G/D614G_all_isolates.csv")


out_file_fasta = "data/9_29_allD614G/D614G_all_isolates_fastaq.csv"

fastas_directory = "data/9_29_allD614G/all_fastas/"
fasta_analysis.getSequenceQuality(fastas_directory,out_file_fasta)

df_fasta = pd.read_csv("data/9_29_allD614G/D614G_all_isolates_fastaq.csv")
df_fasta["Strain"] = [x.replace("-R","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("r1","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("v3","").replace("V3","") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("_","-") for x in df_fasta["Strain"].values]
df_fasta["Strain"] = [x.replace("-0","-") for x in df_fasta["Strain"].values]
df_fasta.drop_duplicates("Strain",inplace=True)


dfjoin = pd.merge(df_vcf,df_fasta,how="outer",on="Strain",indicator=True)

# dfjoin.ix[dfjoin["_merge"]=="left_only",["Nproportion"]]=100.0

dfjoin.to_csv("data/9_29_allD614G/D614G_all_withNpercent.csv",index=False)

####################


#####################
### check vcf table and assign mrn to mcovs for specfic mutations

import pandas as pd

df_table = pd.read_excel("data/10_6/results/S/COVID_SNP_MAP_S_vcf_table_10_6.xlsx")
df_table = df_table[ ( (df_table["Variant status"]=="current only") &  (df_table["Tentative 390"]>1.0))]
query_postions = df_table["Genomic locus"].values

df_vcf = pd.read_excel("data/10_6/results/S/COVID_SNP_MAP_S_vcf_10_6.xlsx")
df_vcf = df_vcf[df_vcf.POS.isin(query_postions)]


df_mrn_mcov = pd.read_csv("reference/4_Curated_MCOV_MRN_Strains.csv")
mrn_mcov = {x:y for x,y in zip(df_mrn_mcov.Strain,df_mrn_mcov.MRN)}
colldate_mcov = {x:y for x,y in zip(df_mrn_mcov.Strain,df_mrn_mcov.COLLECTION_DT)}

df_db = pd.read_excel(ref_dir+"mrn_zipcodes.xlsx",converters={'MRN': str})
dfjoin = pd.merge(df_mrn_mcov,df_db, on="MRN",how="left")
zipcode_mcov = {x:y for x,y in zip(dfjoin.Strain,dfjoin.ZIP)}

temp =[]
for indx,row in df_vcf.iterrows():
    variants = ast.literal_eval(row.Variants_Info)
    for variant in variants.keys():
        if variant == row.ALT:
            isolates = variants[variant][1]
            for isolate in isolates:
                temp.append([row.POS, row.REF,row.ALT, row.ALT_count, 
                isolate, mrn_mcov[isolate],colldate_mcov[isolate],zipcode_mcov[isolate]])

df_result = pd.DataFrame(temp,columns=["POS","REF","ALT","ALT-COUNT","ISOLATE","MRN","COLLECTION-DATE","ZIP"])
df_result.to_excel("postion_mcov_mrn.xlsx",index=False)


df = pd.read_excel("postion_mcov_mrn.xlsx")
df = df[df.POS=24076]
vcf_mcovs = df.ISOLATE.values

results={}
for vcf_mcov in vcf_mcovs:
    df_strain = pd.read_csv("/home/tmhsxs240/COVID_19/data/10_6/10_6_vcfs/"+vcf_mcov+".primertrimmed.medaka.vcf",sep="\t",skiprows=9)
    temp=[]
    for indx,row in df_strain.iterrows():
        temp.append(row.REF+row.POS+row.ALT)
    results[vcf_mcov] = temp


### plot histogram of d614g quality for all isolates

df = pd.read_csv("D614G_all_withmeta.csv")
df = df[df.POS=="23403"]
df = df[df.COLLECTION_DT.notnull()]
df.COLLECTION_DT  = [x[0:10] for x in df.COLLECTION_DT.values]
df['COLLECTION_DT'] = pd.to_datetime(df['COLLECTION_DT'], format='%Y-%m-%d')
df.set_index("COLLECTION_DT",inplace=True)
df_week = df.to_period(freq='W-SAT')
df_week.head()
df_week.reset_index(inplace=True)
df_week.sort_values("COLLECTION_DT",inplace=True)
plt.figure(num=None, figsize=(15, 10), dpi=300, facecolor='w', edgecolor='k')
sns.boxplot(x="COLLECTION_DT", y="QUAL", data=df_week)
plt.xticks(rotation=75)
plt.tight_layout(pad=2.0)
plt.savefig("temp.png");plt.close()
df_week.to_csv("df_week_d614g_qual.csv",index=False)




#### get variants per 500bp for july vs aug-sept
import pandas as pd 

df = pd.read.csv("reference/4_2_Curated_MCOV_MRN_Strains_good_strains.csv")

df_july = df[df.Group=="July"]
df_augsept = df[df.Group=="AugSept"]


df_clade = pd.read_csv("reference/clade_assignment_Oct.csv")

df_july_clade = pd.merge(df_july,df_clade,on="Strain",how="left")
df_july_clade["Clade"].value_counts()

df_augsept_clade = pd.merge(df_augsept,df_clade,on="Strain",how="left")
df_augsept_clade["Clade"].value_counts()

df_gwmut = pd.read_csv("reference/5_StrainsAndMutations_entire_genome.csv")

df_july_gwmut = pd.merge(df_july,df_gwmut,on="Strain",how="left")

df_augsept_gwmut = pd.merge(df_augsept,df_gwmut,on="Strain",how="left")


def fil_position_gwmut(gwmut,df):
    for indx,row in df.iterrows():
        for mutation in ast.literal_eval(row["mutation_info"]):
            mut_position =int(mutation[1:len(mutation)-1])
            gwmut[mut_position] += 1
    return gwmut


july_gwmut={x:0 for x in range(30000)}
augsept_gwmut={x:0 for x in range(30000)}

july_gwmut = fil_position_gwmut(july_gwmut,df_july_gwmut)
augsept_gwmut = fil_position_gwmut(augsept_gwmut,df_augsept_gwmut)

df_july_gwmut_count = pd.DataFrame.from_dict(july_gwmut,orient='index')
df_july_gwmut_count.columns=["count"]
df_july_gwmut_count["count_n"] = df_july_gwmut_count["count"]/df_july_gwmut_count.shape[0]

df_augsept_gwmut_count = pd.DataFrame.from_dict(augsept_gwmut,orient='index')
df_augsept_gwmut_count.columns=["count"]
df_augsept_gwmut_count["count_n"] = df_augsept_gwmut_count["count"]/df_augsept_gwmut_count.shape[0]


july_count = df_july_gwmut_count.groupby(df_july_gwmut_count.index //500).sum()
augsept_count = df_augsept_gwmut_count.groupby(df_augsept_gwmut_count.index //500).sum()


df_july_gwmut_count["count"].plot(color='red')
plt.savefig("july.png");plt.close()

df_augsept_gwmut_count["count"].plot()
plt.savefig("augsept.png");plt.close()


df_july_gwmut_count["count_n"].plot(color='red')
plt.savefig("july.png");plt.close()

df_augsept_gwmut_count["count_n"].plot()
plt.savefig("augsept.png");plt.close()

fig, axes = plt.subplots(2,figsize=(20, 10),dpi=400)
fig.suptitle('Variant Distribution across SARS-CoV-2 Genome')
df_july_gwmut_count["count_n"].plot(color='red',ax=axes[0])
df_augsept_gwmut_count["count_n"].plot(ax=axes[1])
plt.savefig("july_augsept.png");plt.close()


###add region
df_july_gwmut_count["region"]=""
for indx,row in df_july_gwmut_count.iterrows():
    if indx>=266 and indx<=21555:
        df_july_gwmut_count.ix[indx,"region"] = "orf1ab"
    elif indx>=21563 and indx<=25384:
        df_july_gwmut_count.ix[indx,"region"] = "S"
    elif indx>=25393 and indx<=26220:
        df_july_gwmut_count.ix[indx,"region"] = "ORF3a"
    elif indx>=26245 and indx<=26472:
        df_july_gwmut_count.ix[indx,"region"] = "E"
    elif indx>=26523 and indx<=27191:
        df_july_gwmut_count.ix[indx,"region"] = "M"
    elif indx>=27202 and indx<=27387:
        df_july_gwmut_count.ix[indx,"region"] = "ORF6"
    elif indx>=27394 and indx<=27759:
        df_july_gwmut_count.ix[indx,"region"] = "ORF7a"
    elif indx>=27894 and indx<=28259:
        df_july_gwmut_count.ix[indx,"region"] = "ORF8"
    elif indx>=28274 and indx<=29533:
        df_july_gwmut_count.ix[indx,"region"] = "N"
    elif indx>=29558 and indx<=29674:
        df_july_gwmut_count.ix[indx,"region"] = "ORF10"

df_july_gwmut_count["augsept_n"]=df_augsept_gwmut_count["count_n"]
df_july_gwmut_count["augsept"]=df_augsept_gwmut_count["count"]
df_july_gwmut_count.to_csv("july_augsept_data.csv")

